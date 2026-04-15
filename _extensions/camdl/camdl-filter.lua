--- camdl-filter.lua
--- Quarto filter for {camdl} code blocks.
---
--- Extracts .camdl source to the filesystem (when #| filename: is set),
--- runs `camdl check` to validate, and renders as a highlighted code block.

local camdl_bin = os.getenv("CAMDL") or "camdl"

--- Parse cell-level options from #| comments at the top of the block.
local function parse_options(code)
  local opts = {}
  local body_lines = {}
  local in_header = true
  -- Split preserving empty lines (gmatch("[^\r\n]+") swallows them)
  local pos = 1
  while pos <= #code do
    local line_end = code:find("\n", pos) or (#code + 1)
    local line = code:sub(pos, line_end - 1):gsub("\r$", "")
    if in_header and line:match("^#|") then
      local key, val = line:match("^#|%s*(%S+)%s*:%s*(.+)$")
      if key and val then
        opts[key] = val:match("^%s*(.-)%s*$")  -- trim
      end
    else
      in_header = false
      table.insert(body_lines, line)
    end
    pos = line_end + 1
  end
  return opts, table.concat(body_lines, "\n")
end

--- Write file content, creating parent directories as needed.
local function write_file(path, content)
  local dir = path:match("(.+)/[^/]+$")
  if dir then
    os.execute("mkdir -p " .. dir)
  end
  local f = io.open(path, "w")
  if f then
    f:write(content)
    f:write("\n")
    f:close()
    return true
  end
  return false
end

--- Run camdl check on the given file path. Returns (ok, message).
local function check_model(path)
  local cmd = camdl_bin .. " check " .. path .. " 2>&1"
  local handle = io.popen(cmd)
  if not handle then
    return false, "failed to run camdl check"
  end
  local output = handle:read("*a")
  local ok, _, code = handle:close()
  if code == 0 then
    return true, output
  else
    return false, output
  end
end

function CodeBlock(block)
  -- Only process blocks with class "camdl"
  if not block.classes:includes("camdl") then
    return nil
  end

  local opts, body = parse_options(block.text)
  local filename = opts["filename"]

  -- Write to file if filename is specified
  if filename then
    local ok = write_file(filename, body)
    if ok then
      -- Check the model compiles
      local check_ok, check_msg = check_model(filename)
      if not check_ok then
        io.stderr:write("\n[camdl-filter] ERROR compiling " .. filename .. ":\n")
        io.stderr:write(check_msg .. "\n")
        -- Still render but add an error callout
        local err_block = pandoc.Div({
          pandoc.Para({pandoc.Strong({pandoc.Str("camdl check failed:")})}),
          pandoc.CodeBlock(check_msg),
        })
        err_block.classes = {"callout-warning"}

        -- Render original code + error
        if filename then
          block.attributes["filename"] = filename
        end
        block.text = body
        return {block, err_block}
      end
    end
  end

  -- Strip #| options from text, keep classes for syntax highlighting
  if filename then
    block.attributes["filename"] = filename
  end
  block.text = body
  return block
end
