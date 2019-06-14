var calculate;

calculate = function() {
  var ATOMICMASSUNIT, a0, g, re, result, target, task, type, z0;
  ATOMICMASSUNIT = 931.4943;
  task = $('input[name=task]:checked').val();
  re = Number($('#RE').val());
  z0 = Number($('#Z').val());
  a0 = Number($('#A').val());
  target = $('#select_target').val();
  switch (task) {
    case 'r':
      type = 'Range';
      result = 1.23;
      break;
    case 'e':
      type = 'Energy';
      result = 950.333;
      break;
    case 'd':
      g = 1.0 + re / ATOMICMASSUNIT;
      type = 'dE/dx';
      result = g;
  }
  $('#result').html(type + ": " + result);
  return true;
};
