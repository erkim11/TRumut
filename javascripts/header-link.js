/* Make the "TRumut" site title in the header a link to the Introduction page.
   Material renders the title as plain text; the hidden logo anchor already
   carries the correct (relative) home URL, so we reuse it. */
(function () {
  function linkTitle() {
    var header = document.querySelector(".md-header");
    if (!header) return;
    var logo = header.querySelector(".md-header__button.md-logo");
    var topic = header.querySelector(".md-header__topic:first-child .md-ellipsis");
    if (!logo || !topic || topic.closest("a")) return;
    var a = document.createElement("a");
    a.href = logo.getAttribute("href");
    a.title = "Go to Introduction";
    a.className = "md-header__home";
    while (topic.firstChild) a.appendChild(topic.firstChild);
    topic.appendChild(a);
  }
  if (typeof document$ !== "undefined") {
    document$.subscribe(linkTitle);   /* re-run on instant navigation */
  } else {
    document.addEventListener("DOMContentLoaded", linkTitle);
  }
})();
