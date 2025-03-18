// includes short js script to get the correct window size for pdf plotting

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

// Listen for the custom message from the server
Shiny.addCustomMessageHandler('update-sidebar-height', function(message) {
    // Find the plot container by its ID
    var sidebar = document.querySelector(".sidebar"); // Select the sidebar element
    if (sidebar) {
        // Dynamically adjust the sidebar's max-height to match the plot height
        sidebar.style.setProperty("max-height", message.height, "important");
    }
});