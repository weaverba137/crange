/// <reference path="../../../../local/products/node_modules/lib/node_modules/@types/jquery/index.d.ts" />
interface Target {
    name: string;
    z2: number;
    a2: number;
    iadj: number;
    rho: number;
    pla: number;
    etad: number;
    bind: number;
    X0: number;
    X1: number;
    a: number;
    m: number;
    d0: number;
}

$(
    function(): void {
        let targets: Target[] = [];
        let absorberTable = function(): void {
            const ref: string[] = ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0'];
            let body: JQuery<HTMLElement> = $('#atbody');
            for (let i: number = 0; i < targets.length; i++) {
                if (targets[i].name !== 'Unknown') {
                    let rowid: string = 't' + i;
                    $("<tr id=\"" + rowid + "\"/>").appendTo(body);
                    let rowref: JQuery<HTMLElement> = $('#' + rowid);
                    for (let j: number = 0; j < ref.length; j++) {
                        $('<td/>').html(targets[ref[j]]).appendTo(rowref);
                    }
                }
            }
        };
        if (targets.length === 0) {
            $.getJSON("target.json", {}, function(data): void {
                targets = data;
            }).fail(function(): void {
                alert("JSON error!");
            }).done(absorberTable);
        }
    }
);
