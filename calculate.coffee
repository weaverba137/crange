calculate = ->
    ATOMICMASSUNIT = 931.4943
    task = $('input[name=task]:checked').val()
    re = Number $('#RE').val()
    z0 = Number $('#Z').val()
    a0 = Number $('#A').val()
    target = $('#select_target').val()
    switch task
        when 'r'
            type = 'Range'
            result = 1.23
        when 'e'
            type = 'Energy'
            result = 950.333
        when 'd'
            g = 1.0 + re/ATOMICMASSUNIT
            type = 'dE/dx'
            result = g
    $('#result').html("#{type}: #{result}")
    true
