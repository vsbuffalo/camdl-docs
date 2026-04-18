// Add copy buttons to pre.ansi-output blocks (CLI output from run_cli).
// Quarto's code-copy only covers its own rendered code blocks.
document.addEventListener("DOMContentLoaded", function () {
  document.querySelectorAll("pre.ansi-output").forEach(function (pre) {
    var btn = document.createElement("button");
    btn.className = "code-copy-button";
    btn.title = "Copy to clipboard";
    btn.innerHTML =
      '<i class="bi bi-clipboard"></i>';
    btn.addEventListener("click", function () {
      var code = pre.querySelector("code");
      // Strip the leading "$ command" prompt line and HTML tags for plain text
      var text = code.innerText.replace(/^\$ .+\n/, "");
      navigator.clipboard.writeText(text).then(function () {
        btn.innerHTML = '<i class="bi bi-check2"></i>';
        setTimeout(function () {
          btn.innerHTML = '<i class="bi bi-clipboard"></i>';
        }, 1500);
      });
    });
    pre.style.position = "relative";
    btn.style.position = "absolute";
    btn.style.top = "0.3em";
    btn.style.right = "0.3em";
    pre.appendChild(btn);
  });
});
