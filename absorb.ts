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
            const col: string[] = ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0'];
            for (let i: number = 0; i < targets.length; i++) {
                if (targets[i].name !== 'Unknown') {
                    $("<tr id=\"t" + i + "\"/>").appendTo($('#atbody'));
                    for (let j: number = 0; j < col.length; j++) {
                        $('<td/>').html(targets[i][col[j]]).appendTo($('#t' + i));
                    }
                }
            }
        };
        if (targets.length === 0) {
            $.getJSON("target.json", {}, function(data: Target[]): void {
                targets = data;
            }).fail(function(): void {
                alert("JSON error!");
            }).done(absorberTable);
        }
    }
);
