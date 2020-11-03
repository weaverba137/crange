$(function () {
    let targets = [];
    let absorberTable = function () {
        const ref = ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0'];
        let body = $('#atbody');
        for (let i = 0; i < targets.length; i++) {
            if (targets[i].name !== 'Unknown') {
                let rowid = 't' + i;
                $("<tr id=\"" + rowid + "\"/>").appendTo(body);
                let rowref = $('#' + rowid);
                for (let j = 0; j < ref.length; j++) {
                    $('<td/>').html(targets[ref[j]]).appendTo(rowref);
                }
            }
        }
    };
    if (targets.length === 0) {
        $.getJSON("target.json", {}, function (data) {
            targets = data;
        }).fail(function () {
            alert("JSON error!");
        }).done(absorberTable);
    }
});
