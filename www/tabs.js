$(document).ready(function () {
  // initialise state with all tabs hidden
  $("#tabLoad").hide();
  $("#tabAnalysis").hide();
  $("#tabPairwise").hide();
  $("#tabRob").hide();
  $('#splash').hide();
  $('#splash').show();

  // TODO: handle all tabs
  $('#nav-home').click(function() {
    Shiny.unbindAll();
    $('#splash').show();
    $('#tabLoad').hide();
    $('#tabAnalysis').hide();
    $('#tabPairwise').hide();
    $('#tabRob').hide();
    // TODO: remove from all .nav-links
    $('#nav-load').removeClass('active')
    $('#nav-analysis').removeClass('active')
    $('#nav-home').addClass('active');
    $('#nav-pairwise').removeClass('active');
    $('#nav-rob').removeClass('active');
    Shiny.bindAll();
  });

  $('#nav-load').click(function() {
    Shiny.unbindAll();
    $('#splash').hide();
    $('#tabAnalysis').hide();
    $('#tabLoad').show();
    $('#tabPairwise').hide();
    $('#tabRob').hide();
    $('#nav-load').className = 'nav-link';
    // TODO: remove from all .nav-links
    $('#nav-home').removeClass('active')
    $('#nav-analysis').removeClass('active')
    $('#nav-load').addClass('active');
    $('#nav-pairwise').removeClass('active');
    $('#nav-rob').removeClass('active');
    Shiny.bindAll();
  });

  $('#nav-analysis').click(function() {
    Shiny.unbindAll();
    $('#splash').hide();
    $('#tabLoad').hide();
    $('#tabAnalysis').show();
    $('#tabPairwise').hide();
    $('#tabRob').hide();
    // TODO: remove from all .nav-links
    $('#nav-home').removeClass('active')
    $('#nav-load').removeClass('active')
    $('#nav-analysis').addClass('active');
    $('#nav-pairwise').removeClass('active');
    $('#nav-rob').removeClass('active');
    Shiny.bindAll();

    new $.fn.tooltip.Constructor($('html'), {
      "title": "<b>Unrelated</b>: The 'basic' interaction terms are given unrelated," +
               " vague priors<br><b>Exchangeable</b>: The 'basic' interaction terms" +
               " are drawn from a random distribution with a common mean and between-treatment variance",
      "selector": "#inputBeta",
      "placement": "right",
      "html": true,
    });
  });

  $('#nav-pairwise').click(function() {
    Shiny.unbindAll();
    $('#splash').hide();
    $('#tabLoad').hide();
    $('#tabAnalysis').hide();
    $('#tabPairwise').show();
    $('#tabRob').hide();
    // TODO: remove from all .nav-links
    $('#nav-home').removeClass('active')
    $('#nav-load').removeClass('active')
    $('#nav-analysis').removeClass('active');
    $('#nav-pairwise').addClass('active');
    $('#nav-rob').removeClass('active');
    Shiny.bindAll();
  });

  $('#nav-rob').click(function() {
    Shiny.unbindAll();
    $('#splash').hide();
    $('#tabLoad').hide();
    $('#tabAnalysis').hide();
    $('#tabPairwise').hide();
    $('#tabRob').show();
    // TODO: remove from all .nav-links
    $('#nav-home').removeClass('active')
    $('#nav-load').removeClass('active')
    $('#nav-analysis').removeClass('active');
    $('#nav-pairwise').removeClass('active');
    $('#nav-rob').addClass('active');
    Shiny.bindAll();
  });
});

function startAnalysis () {
  var r = confirm("Start Analysis?");
  if (r === true) {
    console.log("starting analysis")
    Shiny.setInputValue("startAnalysis", true);
  } else {
    console.log("Analysis cancelled")
  }
}

function selectPairwiseTableItem(id) {
  var elem = document.getElementById(id).value;
  var out = {id:id, value:elem};
  Shiny.setInputValue("pairwiseTableSelector", out);
}

