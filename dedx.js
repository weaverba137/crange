$(function() {
  var ALPHA, ATOMICMASSUNIT, DeDx, ELECTRONMASS, PROTONMASS, absorber, absorberTable, calculate, dedx, div, effective_charge, i, len, ref, validateNumber;
  DeDx = {
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
    ]
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
  absorber = function(target) {
    var i, len, ref, t;
    ref = DeDx.targets;
    for (i = 0, len = ref.length; i < len; i++) {
      t = ref[i];
      if (t.name === target) {
        return t;
      }
    }
    return t;
  };
  ALPHA = 7.29735301383e-3;
  ATOMICMASSUNIT = 931.4943;
  PROTONMASS = 938.2723;
  ELECTRONMASS = 0.511003e+6;
  effective_charge = function(z0, e1, z2) {
    var b, b2, capA, capB, g, z23;
    g = 1.0 + e1 / ATOMICMASSUNIT;
    b2 = 1.0 - 1.0 / (g * g);
    b = math.sqrt(b2);
    z23 = math.exp((2.0 / 3.0) * math.log(z2));
    capA = 1.16 - z2 * (1.91e-03 - 1.26e-05 * z2);
    capB = (1.18 - z2 * (7.5e-03 - 4.53e-05 * z2)) / ALPHA;
    return z0 * (1.0 - capA * math.exp(-capB * b / z23));
  };
  dedx = function(e1, z0, a1, t) {
    var b, b2, f1, f2, f6, g, z1;
    g = 1.0 + e1 / ATOMICMASSUNIT;
    b2 = 1.0 - 1.0 / (g * g);
    b = math.sqrt(b2);
    z1 = effective_charge(z0, e1, t.z2);
    f1 = 0.3070722 * z1 * z1 * t.z2 / (b2 * a1 * t.a2);
    f2 = math.log(2.0 * ELECTRONMASS * b2 / t.iadj);
    f6 = 2.0 * math.log(g) - b2;
    return f1 * (f2 + f6);
  };
  calculate = function(eventObject) {
    var a1, re, result, target, task, type, unit, z0;
    task = $('input[name=task]:checked').val();
    re = Number($('#RE').val());
    z0 = Number($('#Z').val());
    a1 = Number($('#A').val());
    target = absorber($('#select_target').val());
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
        type = 'dE/dx';
        result = dedx(re, z0, a1, target);
        unit = 'A&nbsp;MeV&nbsp;g<sup>-1</sup>&nbsp;cm<sup>2</sup>';
    }
    $('#result').html(type + ": " + result + " " + unit);
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
  return true;
});
