--- ansi.lua
--- Quarto filter that converts ANSI SGR escape sequences in code blocks
--- to styled HTML <span> elements. Use class "ansi" on fenced code blocks:
---
---   ```{.ansi}
---   error[E303]: conflicting dimensions
---   ```
---
--- Or capture real CLI output by running a command in a {bash} block and
--- piping through this filter.

-- SGR code → CSS style mapping
local sgr_styles = {
  -- Reset
  ["0"]  = nil,
  -- Bold / faint
  ["1"]  = "font-weight:bold",
  ["2"]  = "opacity:0.6",
  -- Italic / underline
  ["3"]  = "font-style:italic",
  ["4"]  = "text-decoration:underline",
  -- Foreground colors (standard)
  ["30"] = "color:#1a1a1a",  -- black
  ["31"] = "color:#cc3333",  -- red
  ["32"] = "color:#2e8b2e",  -- green
  ["33"] = "color:#b8860b",  -- yellow/dark yellow
  ["34"] = "color:#3366cc",  -- blue
  ["35"] = "color:#9b59b6",  -- magenta
  ["36"] = "color:#1abc9c",  -- cyan
  ["37"] = "color:#cccccc",  -- white
  -- Bright foreground
  ["90"] = "color:#666666",  -- bright black (gray)
  ["91"] = "color:#e74c3c",  -- bright red
  ["92"] = "color:#27ae60",  -- bright green
  ["93"] = "color:#f39c12",  -- bright yellow
  ["94"] = "color:#5dade2",  -- bright blue
  ["95"] = "color:#c39bd3",  -- bright magenta
  ["96"] = "color:#48c9b0",  -- bright cyan
  ["97"] = "color:#ffffff",  -- bright white
}

--- Parse a semicolon-separated SGR parameter string into a combined CSS style.
local function sgr_to_css(params)
  local styles = {}
  for code in params:gmatch("[^;]+") do
    local s = sgr_styles[code]
    if s then
      table.insert(styles, s)
    end
  end
  if #styles == 0 then
    return nil
  end
  return table.concat(styles, ";")
end

--- Convert a string containing ANSI escape sequences to HTML with inline styles.
local function ansi_to_html(text)
  local result = {}
  local pos = 1
  local open_spans = 0

  while pos <= #text do
    -- Match ESC[ ... m  (CSI SGR sequence)
    local esc_start, esc_end, params = text:find("\027%[([%d;]*)m", pos)

    if not esc_start then
      -- No more escapes — emit the rest as-is
      table.insert(result, (text:sub(pos):gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;")))
      break
    end

    -- Emit text before this escape
    if esc_start > pos then
      local chunk = text:sub(pos, esc_start - 1)
      table.insert(result, (chunk:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;")))
    end

    -- Process the SGR
    if params == "" or params == "0" then
      -- Reset: close all open spans
      for _ = 1, open_spans do
        table.insert(result, "</span>")
      end
      open_spans = 0
    else
      local css = sgr_to_css(params)
      if css then
        table.insert(result, '<span style="' .. css .. '">')
        open_spans = open_spans + 1
      end
    end

    pos = esc_end + 1
  end

  -- Close any remaining open spans
  for _ = 1, open_spans do
    table.insert(result, "</span>")
  end

  return table.concat(result)
end

function CodeBlock(block)
  if not block.classes:includes("ansi") then
    return nil
  end

  -- Only produce HTML output
  if not quarto.doc.is_format("html:js") then
    -- For non-HTML formats, strip ANSI codes and return plain text
    local plain = block.text:gsub("\027%[%d*;?%d*m", "")
    block.text = plain
    block.classes = {"text"}
    return block
  end

  local html = ansi_to_html(block.text)

  -- Wrap in a <pre><code> block styled like a code block
  local wrapped = '<pre class="ansi-output"><code>' .. html .. '</code></pre>'

  return pandoc.RawBlock("html", wrapped)
end
