$(function() {
  var DeDx, absorberTable, calculate, csv, div, i, len, ref, resetForm, validateNumber;
  DeDx = {
    nCalc: 0,
    ALPHA: 7.29735301383e-3,
    ATOMICMASSUNIT: 931.4943,
    PROTONMASS: 938.2723,
    ELECTRONMASS: 0.511003e+6,
    targets: [],
    validateDivs: [
      {
        name: 'div_task',
        type: 'radio',
        valid: false,
        first: true,
        help: ''
      }, {
        name: 'div_re',
        type: 'float',
        valid: false,
        first: true,
        help: 'Units are A&nbsp;MeV or g&nbsp;cm<sup>-2</sup>.'
      }, {
        name: 'div_z',
        type: 'int',
        valid: false,
        first: true,
        help: ''
      }, {
        name: 'div_a',
        type: 'int',
        valid: false,
        first: true,
        help: ''
      }
    ],
    absorber: function(target) {
      var i, len, ref, t;
      ref = this.targets;
      for (i = 0, len = ref.length; i < len; i++) {
        t = ref[i];
        if (t.name === target) {
          return t;
        }
      }
      return t;
    },
    effective_charge: function(z0, e1, z2) {
      var b, b2, capA, capB, g, z23;
      g = 1.0 + e1 / this.ATOMICMASSUNIT;
      b2 = 1.0 - 1.0 / (g * g);
      b = math.sqrt(b2);
      z23 = math.exp((2.0 / 3.0) * math.log(z2));
      capA = 1.16 - z2 * (1.91e-03 - 1.26e-05 * z2);
      capB = (1.18 - z2 * (7.5e-03 - 4.53e-05 * z2)) / this.ALPHA;
      return z0 * (1.0 - capA * math.exp(-capB * b / z23));
    },
    delta: function(g, t) {
      var X, X0, X1, b, cbar;
      X0 = t.X0;
      X1 = t.X1;
      cbar = 2.0 * math.log(t.iadj / t.pla) + 1.0;
      b = math.sqrt(1.0 - 1.0 / (g * g));
      X = math.log10(b * g);
      if (t.etad > 0) {
        cbar -= 2.303 * math.log10(t.etad);
        X1 -= 0.5 * math.log10(t.etad);
        X0 -= 0.5 * math.log10(t.etad);
      }
      if (X < X0) {
        return t.d0 * math.exp(4.6052 * (X - X0));
      }
      if (X >= X0 && X < X1) {
        return 4.6052 * X + math.exp(math.log(t.a) + t.m * math.log(X1 - X)) - cbar;
      }
      return 4.6052 * X - cbar;
    },
    dedx: function(e1, z0, a1, target) {
      var b, b2, delt, f1, f2, f6, g, t, z1;
      t = this.absorber(target);
      g = 1.0 + e1 / this.ATOMICMASSUNIT;
      delt = this.delta(g, t);
      b2 = 1.0 - 1.0 / (g * g);
      b = math.sqrt(b2);
      z1 = this.effective_charge(z0, e1, t.z2);
      f1 = 0.3070722 * z1 * z1 * t.z2 / (b2 * a1 * t.a2);
      f2 = math.log(2.0 * this.ELECTRONMASS * b2 / t.iadj);
      f6 = 2.0 * math.log(g) - b2;
      return f1 * (f2 + f6 + (delt / 2.0));
    }
  };
  validateNumber = function(eventObject) {
    var d, div, formvalid, input, name, patt, previous_form, validity;
    previous_form = $('#previous_target').val() !== 'Unknown';
    if (previous_form) {
      eventObject.data.first = false;
    }
    div = eventObject.data.name;
    if (eventObject.data.type === 'radio') {
      name = div.split('_')[1];
      eventObject.data.valid = $("input[name=" + name + "]:checked").length === 1;
      if (!eventObject.data.first) {
        if (eventObject.data.valid) {
          $('#' + div).removeClass('has-error');
        } else {
          $('#' + div).addClass('has-error');
        }
      }
    } else {
      input = $('#' + div.split('_')[1].toUpperCase());
      patt = eventObject.data.type === 'int' ? new RegExp(/^[0-9]+$/i) : new RegExp(/^[0-9]+(\.[0-9]*|)(e[+-]?[0-9]+|)$/i);
      eventObject.data.valid = patt.test(input.val());
      if (!eventObject.data.first) {
        if (eventObject.data.valid) {
          $('#' + div).removeClass('has-error').addClass('has-success');
          $('#' + div + '_helpblock').html(eventObject.data.help);
        } else {
          $('#' + div).removeClass('has-success').addClass('has-error');
          $('#' + div + '_helpblock').html(eventObject.data.help + " " + (eventObject.data.type === 'int' ? 'Integer' : 'Float') + " value required!");
        }
      }
    }
    eventObject.data.first = false;
    validity = (function() {
      var i, len, ref, results;
      ref = DeDx.validateDivs;
      results = [];
      for (i = 0, len = ref.length; i < len; i++) {
        d = ref[i];
        results.push(d.valid);
      }
      return results;
    })();
    formvalid = validity.every(function(currentValue) {
      return currentValue;
    });
    $('#calculate').prop('disabled', !formvalid);
    return eventObject.data.valid;
  };
  absorberTable = function() {
    var i, k, len, previous_target, ref, rowid, select, selected_target, t;
    select = $('#select_target');
    previous_target = $('#previous_target');
    selected_target = previous_target.length === 1 ? previous_target.val() : 'Unknown';
    k = 0;
    ref = DeDx.targets;
    for (i = 0, len = ref.length; i < len; i++) {
      t = ref[i];
      if (t.name !== 'Unknown') {
        rowid = 't' + k;
        k += 1;
        $("<option id=\"" + rowid + "\"/>").html(t.name).appendTo(select);
        if (t.name === selected_target) {
          $("#" + rowid).prop('selected', true);
        }
      }
    }
    return k;
  };
  calculate = function(eventObject) {
    var a1, r, re, result, rr, target, task, type, unit, z0;
    task = $('input[name=task]:checked').val();
    re = Number($('#RE').val());
    z0 = Number($('#Z').val());
    a1 = Number($('#A').val());
    target = $('#select_target').val();
    switch (task) {
      case 'r':
        type = 'Range';
        result = 1.23;
        unit = 'g&nbsp;cm<sup>-2</sup>';
        break;
      case 'e':
        type = 'Energy';
        result = 950.333;
        unit = 'A&nbsp;MeV';
        break;
      case 'd':
        type = 'd<var>E</var>/d<var>x</var>';
        result = DeDx.dedx(re, z0, a1, target);
        unit = 'A&nbsp;MeV&nbsp;g<sup>-1</sup>&nbsp;cm<sup>2</sup>';
    }
    DeDx.nCalc += 1;
    r = $('#result');
    rr = $("<tr id=\"rtr" + DeDx.nCalc + "\"/>");
    $('<td/>').html(DeDx.nCalc).appendTo(rr);
    $('<td/>').html(type).appendTo(rr);
    $('<td/>').html(re).appendTo(rr);
    $('<td/>').html(z0).appendTo(rr);
    $('<td/>').html(a1).appendTo(rr);
    $('<td/>').html(target).appendTo(rr);
    $('<td/>').html(result).appendTo(rr);
    $('<td/>').html(unit).appendTo(rr);
    rr.appendTo(r);
    return true;
  };
  csv = function(eventObject) {
    var c, col, download, foo, header, i, j, len, len1, r, ref, row, rows;
    rows = [];
    header = ['ID', 'Task', 'E/R', 'Z', 'A', 'Target', 'Result', 'Units'];
    rows.push(header.join(','));
    foo = $('#result').children();
    if (foo.length === 0) {
      alert('No rows!');
      return false;
    }
    for (i = 0, len = foo.length; i < len; i++) {
      row = foo[i];
      r = [];
      ref = row.children;
      for (j = 0, len1 = ref.length; j < len1; j++) {
        col = ref[j];
        r.push(col.innerHTML.replace(/&nbsp;/g, ' '));
      }
      rows.push(r.join(','));
    }
    c = rows.join('\r\n') + '\r\n';
    download = $("<a/>", {
      href: 'data:text/csv;charset=utf-8,' + encodeURIComponent(c),
      download: "crange.csv"
    });
    download.appendTo($("#downloadCSV"));
    download[0].click();
    $("#downloadCSV").empty();
    return true;
  };
  resetForm = function(eventObject) {
    document.getElementById('recalc').reset();
    $('#result').empty();
    DeDx.nCalc = 0;
    return true;
  };
  if (DeDx.targets.length === 0) {
    $.getJSON("target.json", {}, function(data) {
      return DeDx.targets = data;
    }).fail(function() {
      return alert("JSON error!");
    }).done(absorberTable);
  }
  ref = DeDx.validateDivs;
  for (i = 0, len = ref.length; i < len; i++) {
    div = ref[i];
    if ($("#" + div.name).length > 0) {
      $("#" + div.name).change(div, validateNumber).change();
    }
  }
  $('#calculate').click(calculate);
  $('#resetForm').click(resetForm);
  $('#CSV').click(csv);
  return true;
});
