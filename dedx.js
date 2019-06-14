$(function() {
  var ALPHA, ATOMICMASSUNIT, DeDx, ELECTRONMASS, PROTONMASS, absorber, absorberTable, calculate, dedx, div, i, len, ref, validateNumber;
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
    $('#engage').prop('disabled', !formvalid);
    return eventObject.data.valid;
  };
  absorberTable = function() {
    var body, i, j, k, l, len, len1, previous_target, ref, ref1, rowid, rowref, select, selected_target, t;
    body = $('#atbody');
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
        if (body.length === 1) {
          $("<tr id=\"" + rowid + "\"/>").appendTo(body);
          rowref = $('#' + rowid);
          ref1 = ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0'];
          for (j = 0, len1 = ref1.length; j < len1; j++) {
            l = ref1[j];
            $('<td/>').html(t[l]).appendTo(rowref);
          }
        }
        if (select.length === 1) {
          $("<option id=\"" + rowid + "\"/>").html(t.name).appendTo(select);
          if (t.name === selected_target) {
            $("#" + rowid).prop('selected', true);
          }
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
  dedx = function(e1, z0, a1, t) {
    var b, b2, f1, f2, f6, g, z1;
    g = 1.0 + e1 / ATOMICMASSUNIT;
    b2 = 1.0 - 1.0 / (g * g);
    b = math.sqrt(b2);
    z1 = z0;
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
  $('#engage').click(calculate);
  return true;
});
