$ () ->
    DeDx =
        targets: []
        validateDivs: [
            {name: 'div_task', type: 'radio', valid: false, first: true, help: ''}
            {name: 'div_re', type: 'float', valid: false, first: true, help: 'Units are A&nbsp;MeV or g&nbsp;cm<sup>-2</sup>.'}
            {name: 'div_z', type: 'int', valid: false, first: true, help: ''}
            {name: 'div_a', type: 'int', valid: false, first: true, help: ''}
            ]
    #
    # Validate numerical values
    #
    validateNumber = (eventObject) ->
        previous_form = $('#previous_target').val() != 'Unknown'
        eventObject.data.first = false if previous_form
        div = eventObject.data.name
        if eventObject.data.type == 'radio'
            name = div.split('_')[1]
            eventObject.data.valid = $("input[name=#{name}]:checked").length == 1
            if not eventObject.data.first
                if eventObject.data.valid
                    $('#'+div).removeClass('has-error')
                else
                    $('#'+div).addClass('has-error')

        else
            input = $('#' + div.split('_')[1].toUpperCase())
            patt = if eventObject.data.type == 'int' then new RegExp /^[0-9]+$/i else new RegExp /^[0-9]+(\.[0-9]*|)(e[+-]?[0-9]+|)$/i
            eventObject.data.valid = patt.test input.val()
            if not eventObject.data.first
                if eventObject.data.valid
                    $('#'+div).removeClass('has-error').addClass('has-success')
                    $('#'+div+'_helpblock').html(eventObject.data.help)
                else
                    $('#'+div).removeClass('has-success').addClass('has-error')
                    $('#'+div+'_helpblock').html("#{eventObject.data.help} #{if eventObject.data.type == 'int' then 'Integer' else 'Float'} value required!")
        eventObject.data.first = false
        validity = (d.valid for d in DeDx.validateDivs)
        formvalid = validity.every (currentValue) -> currentValue
        $('#calculate').prop('disabled', not formvalid)
        eventObject.data.valid
    #
    # Convert target data into table.
    #
    absorberTable = ->
        select = $('#select_target')
        previous_target = $('#previous_target')
        selected_target = if previous_target.length == 1 then previous_target.val() else 'Unknown'
        k = 0
        for t in DeDx.targets
            if t.name != 'Unknown'
                rowid = 't' + k
                k += 1
                $("<option id=\"#{rowid}\"/>").html(t.name).appendTo select
                $("##{rowid}").prop('selected',true) if t.name == selected_target
        k
    #
    # Get absorber data from target name.
    #
    absorber = (target) ->
        for t in DeDx.targets
            if t.name == target
                return t
        t
    #
    # dE/dx
    #
    ALPHA = 7.29735301383e-3
    ATOMICMASSUNIT = 931.4943
    PROTONMASS = 938.2723
    ELECTRONMASS = 0.511003e+6
    effective_charge = (z0, e1, z2) ->
        g = 1.0 + e1/ATOMICMASSUNIT
        b2 = 1.0 - 1.0/(g*g)
        b = math.sqrt(b2)
        z23 = math.exp((2.0/3.0)*math.log(z2))
        capA = 1.16 - z2*(1.91e-03 - 1.26e-05*z2)
        capB = (1.18 - z2*(7.5e-03 - 4.53e-05*z2))/ALPHA
        z0 * (1.0 - capA*math.exp(-capB*b/z23))
    #
    #
    #
    dedx = (e1, z0, a1, t) ->
        g = 1.0 + e1/ATOMICMASSUNIT
        b2 = 1.0 - 1.0/(g*g)
        b = math.sqrt(b2)
        z1 = effective_charge(z0, e1, t.z2)
        f1 = 0.3070722*z1*z1*t.z2/(b2*a1*t.a2)
        f2 = math.log(2.0*ELECTRONMASS*b2/t.iadj)
        f6 = 2.0*math.log(g) - b2
        f1*(f2 + f6)
    #
    # Calculate
    #
    calculate = (eventObject) ->
        task = $('input[name=task]:checked').val()
        re = Number $('#RE').val()
        z0 = Number $('#Z').val()
        a1 = Number $('#A').val()
        target = absorber $('#select_target').val()
        switch task
            when 'r'
                type = 'Range'
                result = 1.23
                unit = 'g&nbsp;cm<sup>-2</sup>'
            when 'e'
                type = 'Energy'
                result = 950.333
                unit = 'A&nbsp;MeV'
            when 'd'
                type = 'dE/dx'
                result = dedx(re, z0, a1, target)
                unit = 'A&nbsp;MeV&nbsp;g<sup>-1</sup>&nbsp;cm<sup>2</sup>'
        $('#result').html("#{type}: #{result} #{unit}")
        true
    #
    # Load the target data.
    #
    if DeDx.targets.length == 0
        $.getJSON("target.json",{},(data) -> DeDx.targets = data).fail(
            () -> alert("JSON error!")
            ).done(absorberTable)
    #
    # Initially disable the Calculate button
    #
    # if $('#calculate').length > 0
    #     $('#calculate').prop('disabled',true)
    #     validateNumber {data: DeDx.validateDivs[0]}
    #
    # Bind to various changes.
    #
    for div in DeDx.validateDivs
        if $("##{div.name}").length > 0
            $("##{div.name}").change(div, validateNumber).change()
    $('#calculate').click(calculate)
    true
