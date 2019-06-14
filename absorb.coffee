$ () ->
    targets = []
    #
    # Convert target data into table.
    #
    absorberTable = ->
        body = $('#atbody')
        k = 0
        for t in targets
            if t.name != 'Unknown'
                rowid = 't' + k
                k += 1
                $("<tr id=\"#{rowid}\"/>").appendTo body
                rowref = $('#'+rowid)
                for l in ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0']
                    $('<td/>').html(t[l]).appendTo rowref
        k
    #
    # Load the target data.
    #
    if targets.length == 0
        $.getJSON("target.json",{},(data) -> targets = data).fail(
            () -> alert("JSON error!")
            ).done(absorberTable)
    true
