$(function () {
    $('[data-toggle="tooltip"]').tooltip()
})

$(function(){
    $('.ali-scroll-div').scroll(function(){
        $('.ali-main-div')
            .scrollLeft($('.ali-scroll-div').scrollLeft());
    });
    $('.ali-main-div').scroll(function(){
        $('.ali-scroll-div')
            .scrollLeft($('.ali-main-div').scrollLeft());
    });
});

function ScoreBreakdown (protein_conformation, cutoff) {
    $.ajax({
        'url': '/alignment/score_breakdown',
        'data': {
            protein_conformation: protein_conformation,
            cutoff: cutoff,
        },
        'type': 'GET',
        'async': false,
        'success': function(data) {
            $("#pconf-" + protein_conformation).html(data);
        },
    });
};

function ApplyCutoff (cutoff) {
    var row = [];
    $('#first-row td').each(function(){
        row.push($(this).attr('id'));
    });

    // show all when hidden
    $('#signature-table tr td[style*="display: none"]').each(function(){
        $(this).css("display", "table-cell");
    });

    // hide when not meeting the cutoff
    var hide = [];
    for (i=0; i < row.length; i++) {
        if (row[i] == 'anchor')
            continue;

        // Only hide when previously shown
        var cell = Math.abs(row[i].substr(7));
        if (cell < cutoff)
            hide.push(i);
    }

    // now hide/show
    $('#signature-table tr').each(function(){
      if ($(this)[0].children.length==row.length) {
        for (i=0; i < hide.length; i++) {
            $(this)[0].children[hide[i]].style.display = "none";
        }
      }
    });

    $('#signature-table #segments td:not("#anchor")').each(function(){
        var segment_name = $(this).attr('id').replace('segment-', '').replace(".","\\.");
        var colspan = $('#gns td#gn-'+segment_name+':visible').length;
        $(this).attr('colspan', colspan);
    });
};

$('#cutoff-apply').click( function() {
    var cutoff = parseInt($('#cutoff-val').val());
    ApplyCutoff(cutoff);
});

/*$(window).on("load", function () {
    $('.internal-scroll-div').css('width', $('.dynamic-div').outerWidth() );
    ApplyCutoff(40);
});*/
