// get the correct window size for pdf plotting
$(document).on('shiny:connected', function(e) {
    var dimensions = {
        width: $(window).width(),
        height: $(window).height()
    };
    Shiny.setInputValue("window_dimensions", dimensions);
});

$(window).on('resize', function(e) {
    var dimensions = {
        width: $(window).width(),
        height: $(window).height()
    };
    Shiny.setInputValue("window_dimensions", dimensions);
});


// Listen for the custom message from the server
Shiny.addCustomMessageHandler("update-plot-height", function(message) {
    // Find the plot container by its ID
    var plotContainer = document.getElementById("msa_plot");
    if (plotContainer) {
        // Update the CSS height property for the plot container
        plotContainer.style.setProperty("height", message.height, "important");
    }
});

// Navbar hides after a while
document.addEventListener("DOMContentLoaded", function () {
    let lastScrollY = window.scrollY;
    const navbar = document.querySelector(".navbar");
    let hideTimeout;
    let isNavbarHidden = false;

    window.addEventListener("scroll", function () {
        clearTimeout(hideTimeout);

        if (window.scrollY > lastScrollY && window.scrollY > 50) {
            if (!isNavbarHidden) {
                hideTimeout = setTimeout(() => {
                    navbar.style.transform = "translateY(-100%)";
                    isNavbarHidden = true;
                }, 1000);
            }
        } else if (window.scrollY < lastScrollY && window.scrollY < 100) {
            navbar.style.transform = "translateY(0)";
            isNavbarHidden = false;
        }

        lastScrollY = window.scrollY;
    });
});

// custom sidebar behaviour
function toggleSidebar() {
    var sidebar = document.getElementById("overlay-sidebar");
    var bg = document.getElementById("overlay-bg");
    if (sidebar.classList.contains("show")) {
        sidebar.classList.remove("show");
        bg.style.display = "none";
    } else {
        sidebar.classList.add("show");
        bg.style.display = "block";
    }
}

// collapse navbar on mobile devices
document.addEventListener("click", function(event) {
    var navbar = document.querySelector(".navbar-collapse");
    var toggle = document.querySelector(".navbar-toggler");

    if (!navbar.contains(event.target) && !toggle.contains(event.target)) {
        navbar.classList.remove("show");  // Bootstrap collapse class
    }
});