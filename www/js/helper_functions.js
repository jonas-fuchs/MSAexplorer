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
    // Flag to check if the navbar is already hidden
    let isNavbarHidden = false;

    window.addEventListener("scroll", function () {
        clearTimeout(hideTimeout);
        // If scrolling down passed the threshold (50px)
        if (window.scrollY > lastScrollY && window.scrollY > 50) {
            // Start the timeout to hide navbar
            if (!isNavbarHidden) {
                hideTimeout = setTimeout(() => {
                    navbar.style.transform = "translateY(-100%)"; // Hide the navbar
                    isNavbarHidden = true; // Update the flag to mark it as hidden
                }, 1000);
            }
        }
        // If scrolling up, show the navbar immediately
        else if (window.scrollY < lastScrollY) {
            navbar.style.transform = "translateY(0)"; // Show the navbar immediately
            isNavbarHidden = false; // Update the flag to mark it as visible
        }

        // Update lastScrollY to the current scroll position
        lastScrollY = window.scrollY;
    });
});
