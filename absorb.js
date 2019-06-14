$(function() {
  var absorberTable, targets;
  targets = [];
  absorberTable = function() {
    var body, i, j, k, l, len, len1, ref, rowid, rowref, t;
    body = $('#atbody');
    k = 0;
    for (i = 0, len = targets.length; i < len; i++) {
      t = targets[i];
      if (t.name !== 'Unknown') {
        rowid = 't' + k;
        k += 1;
        $("<tr id=\"" + rowid + "\"/>").appendTo(body);
        rowref = $('#' + rowid);
        ref = ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0'];
        for (j = 0, len1 = ref.length; j < len1; j++) {
          l = ref[j];
          $('<td/>').html(t[l]).appendTo(rowref);
        }
      }
    }
    return k;
  };
  if (targets.length === 0) {
    $.getJSON("target.json", {}, function(data) {
      return targets = data;
    }).fail(function() {
      return alert("JSON error!");
    }).done(absorberTable);
  }
  return true;
});
