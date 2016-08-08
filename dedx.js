$(function() {
  var DeDx, absorberTable, div, validateNumber, _i, _len, _ref;
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
          $('#' + div + '_helpblock').html("" + eventObject.data.help + " " + (eventObject.data.type === 'int' ? 'Integer' : 'Float') + " value required!");
        }
      }
    }
    eventObject.data.first = false;
    validity = (function() {
      var _i, _len, _ref, _results;
      _ref = DeDx.validateDivs;
      _results = [];
      for (_i = 0, _len = _ref.length; _i < _len; _i++) {
        d = _ref[_i];
        _results.push(d.valid);
      }
      return _results;
    })();
    formvalid = validity.every(function(currentValue) {
      return currentValue;
    });
    $('#engage').prop('disabled', !formvalid);
    return eventObject.data.valid;
  };
  absorberTable = function() {
    var body, k, l, previous_target, rowid, rowref, select, selected_target, t, _i, _j, _len, _len1, _ref, _ref1;
    body = $('#atbody');
    select = $('#select_target');
    previous_target = $('#previous_target');
    selected_target = previous_target.length === 1 ? previous_target.val() : 'Unknown';
    k = 0;
    _ref = DeDx.targets;
    for (_i = 0, _len = _ref.length; _i < _len; _i++) {
      t = _ref[_i];
      if (t.name !== 'Unknown') {
        rowid = 't' + k;
        k += 1;
        if (body.length === 1) {
          $("<tr id=\"" + rowid + "\"/>").appendTo(body);
          rowref = $('#' + rowid);
          _ref1 = ['name', 'z2', 'a2', 'iadj', 'rho', 'pla', 'etad', 'bind', 'X1', 'X1', 'a', 'm', 'd0'];
          for (_j = 0, _len1 = _ref1.length; _j < _len1; _j++) {
            l = _ref1[_j];
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
  if (DeDx.targets.length === 0) {
    $.getJSON("target.json", {}, function(data) {
      return DeDx.targets = data;
    }).error(function() {
      return alert("JSON error!");
    }).complete(absorberTable);
  }
  _ref = DeDx.validateDivs;
  for (_i = 0, _len = _ref.length; _i < _len; _i++) {
    div = _ref[_i];
    if ($("#" + div.name).length > 0) {
      $("#" + div.name).change(div, validateNumber).change();
    }
  }
  return true;
});
