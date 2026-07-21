/* ============================================================
   Renke Wang — Personal Academic Website
   script.js
   Minimal JS: mobile nav toggle, year stamp, active-link highlight.
============================================================ */

(function () {
  "use strict";

  // ---------- Mobile navigation toggle ----------
  var toggle = document.querySelector(".nav-toggle");
  var nav = document.querySelector(".primary-nav");

  if (toggle && nav) {
    toggle.addEventListener("click", function () {
      var isOpen = nav.classList.toggle("open");
      toggle.setAttribute("aria-expanded", String(isOpen));
    });

    // Close menu when a nav link is clicked (mobile UX)
    var links = nav.querySelectorAll("a");
    for (var i = 0; i < links.length; i++) {
      links[i].addEventListener("click", function () {
        nav.classList.remove("open");
        toggle.setAttribute("aria-expanded", "false");
      });
    }

    // Close menu on Escape
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape" && nav.classList.contains("open")) {
        nav.classList.remove("open");
        toggle.setAttribute("aria-expanded", "false");
        toggle.focus();
      }
    });
  }

  // ---------- Double-click a demo video to toggle fullscreen ----------
  // The videos render as a small 3-up grid; visitors double-click any
  // one to watch the animation fullscreen (native controls also keep
  // their own fullscreen button bottom-right).
  var demoVideos = document.querySelectorAll(".demo video");
  for (var d = 0; d < demoVideos.length; d++) {
    demoVideos[d].addEventListener("dblclick", function () {
      var v = this;
      var fsEl = document.fullscreenElement || document.webkitFullscreenElement;
      if (fsEl) {
        (document.exitFullscreen || document.webkitExitFullscreen).call(document);
      } else {
        (v.requestFullscreen || v.webkitRequestFullscreen ||
          v.webkitEnterFullscreen).call(v);
      }
    });
  }

  // ---------- Auto-update copyright year ----------
  var yearEl = document.getElementById("year");
  if (yearEl) {
    yearEl.textContent = String(new Date().getFullYear());
  }

  // ---------- Active section highlighting in nav ----------
  var sections = document.querySelectorAll("main section[id]");
  var navLinks = document.querySelectorAll(".primary-nav a[href^='#']");
  if (sections.length && navLinks.length && "IntersectionObserver" in window) {
    var observer = new IntersectionObserver(
      function (entries) {
        entries.forEach(function (entry) {
          if (entry.isIntersecting) {
            var id = entry.target.getAttribute("id");
            for (var i = 0; i < navLinks.length; i++) {
              var href = navLinks[i].getAttribute("href");
              if (href === "#" + id) {
                navLinks[i].classList.add("active");
              } else {
                navLinks[i].classList.remove("active");
              }
            }
          }
        });
      },
      { rootMargin: "-40% 0px -55% 0px", threshold: 0 }
    );

    for (var j = 0; j < sections.length; j++) {
      observer.observe(sections[j]);
    }
  }
})();
