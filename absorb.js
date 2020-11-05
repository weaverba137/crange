$(function () {
    let targets = [];
    let absorberTable = function () {
        const col = ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0'];
        for (let i = 0; i < targets.length; i++) {
            if (targets[i].name !== 'Unknown') {
                $("<tr id=\"t" + i + "\"/>").appendTo($('#atbody'));
                for (let j = 0; j < col.length; j++) {
                    $('<td/>').html(targets[col[j]]).appendTo($('#t' + i));
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
