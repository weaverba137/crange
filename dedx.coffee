$ () ->
    DeDx =
        nCalc: 0
        ALPHA: 7.29735301383e-3
        ATOMICMASSUNIT: 931.4943
        PROTONMASS: 938.2723
        ELECTRONMASS: 0.511003e+6
        fva: [0.33, 0.078, 0.03, 0.014, 0.0084, 0.0053, 0.0035, 0.0025, 0.0019, 0.0014]
        targets: []
        switches:
            Barkas: false
            Shell: false
            Leung: false
            NewDelta: true
            NewElectronCapture: false
            FiniteNuclearSize: true
            Kinematic: false
            Radiative: false
            Pair: false
            Bremsstrahlung: false
        validateDivs: [
            {name: 'div_task', type: 'radio', valid: false, first: true, help: ''}
            {name: 'div_re', type: 'float', valid: false, first: true, help: 'Units are A&nbsp;MeV or g&nbsp;cm<sup>-2</sup>.'}
            {name: 'div_z', type: 'int', valid: false, first: true, help: ''}
            {name: 'div_a', type: 'int', valid: false, first: true, help: ''}
            ]
        #
        # Get absorber data from target name.
        #
        absorber: (target) ->
            for t in @targets
                if t.name == target
                    return t
            t
        #
        # Effective projectile charge.
        #
        effective_charge: (z0, e1, z2) ->
            g = 1.0 + e1/@ATOMICMASSUNIT
            b2 = 1.0 - 1.0/(g*g)
            b = math.sqrt(b2)
            if @switches.NewElectronCapture
                if z2 == 4.0
                   f1 = (2.045 + 2.000*math.exp(-0.04369*z0))*math.exp(-7.000*math.exp(0.2643*math.log(e1))*math.exp(-0.4171*math.log(z0)))
                else
                    if z2 == 6.0
                        f1 = (2.584 + 1.910*math.exp(-0.03958*z0))*math.exp(-6.933*math.exp(0.2433*math.log(e1))*math.exp(-0.3969*math.log(z0)))
                    else
                        f1 = ((1.164 + 0.2319*math.exp(-0.004302*z2)) + 1.658*math.exp(-0.05170*z0))*math.exp(-(8.144+0.9876*math.log(z2))*math.exp((0.3140+0.01072*math.log(z2))*math.log(e1))*math.exp(-(0.5218+0.02521*math.log(z2))*math.log(z0)))

            else
                z23 = math.exp((2.0/3.0)*math.log(z2))
                #
                # Anthony & Landford:
                #
                capA = 1.16 - z2*(1.91e-03 - 1.26e-05*z2)
                capB = (1.18 - z2*(7.5e-03 - 4.53e-05*z2))/@ALPHA
                #
                # Pierce & Blann:
                #
                # capA = 1.0
                # capB = 130.0
                f1 = capA*math.exp(-capB*b/z23)
            (1.0 - f1)*z0
        #
        # Density effect, Sternheimer, Berger \& Seltzer
        #
        delta: (g, t) ->
            X0 = t.X0
            X1 = t.X1
            cbar = 2.0*math.log(t.iadj/t.pla) + 1.0
            b = math.sqrt(1.0 - 1.0/(g*g))
            X = math.log10(b*g)
            if t.etad > 0
                cbar -= 2.303*math.log10(t.etad)
                X1 -= 0.5*math.log10(t.etad)
                X0 -= 0.5*math.log10(t.etad)
            if X < X0
                return t.d0*math.exp(4.6052*(X-X0))
            if X >= X0 and X < X1
                return (4.6052*X + math.exp(math.log(t.a) + t.m*math.log(X1-X)) - cbar)
            4.6052*X - cbar
        #
        # Obsolete density effect, Sternheimer & Peierls
        #
        olddelta: (g, t) ->
            if g < 1.8
                return 0.0
            cbar = 2.0 * math.log(t.iadj/t.pla)+ 1.0
            b = math.sqrt(1.0 - 1.0/(g*g))
            y = 2.0*math.log10(b*g)
            if t.etad > 0
                y += math.log(etad)
                if cbar >= 12.25
                    y1 = 23.03
                    y0 = if  cbar >= 13.804 then 1.502*cbar-11.52 else 9.212
                else
                    y1 = 18.42
                    if cbar < 12.25
                        y0 = 9.212
                    if cbar < 11.5
                        y0 = 8.751
                    if cbar < 11.0
                        y0 = 8.291
                    if cbar < 10.5
                        y0 = 7.830
                    if cbar < 10.0
                        y0 = 7.370
            else
                if t.iadj >= 100.0
                    y1 = 13.82
                    y0 = if cbar >= 5.215 then 1.502*cbar-6.909 else 0.9212
                else
                    y1 = 9.212
                    y0 = if cbar >= 3.681 then 1.502*cbar-4.606 else 0.9212
            if y < y0
                return 0.0
            else
                if y > y1
                    return y - cbar
                else
                    dy3 = (y1 - y0) * (y1 - y0) * (y1 - y0)
                    a = (cbar - y0)/dy3
                    return (y - cbar + a*(y1 - y) * (y1 - y) * (y1 - y))
        #
        # Complex log of complex Gamma function
        #
        lngamma: (z) ->
            coeff = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
            if z.re > 0
                x = z.re - 1.0
                y = z.im
            else
                x = -z.re
                y = -z.im
            r = math.sqrt((x+5.5)*(x+5.5)+y*y)
            aterm1 = y*math.log(r)
            aterm2 = (x+0.5)*math.atan2(y,(x+5.5))-y
            lterm1 = (x+0.5)*math.log(r)
            lterm2 = -y*math.atan2(y,(x+5.5)) - (x+5.5) + 0.5*math.log(2.0*math.pi)
            num=0.0
            denom=1.000000000190015
            for c, i in coeff
                cterm = c/((x+i+1)*(x+i+1) + y*y)
                num += cterm
                denom += (x + i + 1)*cterm
            num *= -y
            aterm3 = math.atan2(num,denom)
            lterm3 = 0.5*math.log(num*num + denom*denom)
            result = math.complex(lterm1+lterm2+lterm3, aterm1+aterm2+aterm3)
            if z.re < 0
                lpi = math.complex(math.log(math.pi), 0)
                result = math.subtract(lpi, math.add(result, math.log(math.sin(math.multiply(math.pi, z)))))
            result
        #
        # Confluent hypergeometric function
        #
        hyperg: (a, b, z) ->
            dm = 0.0
            term = math.complex(1.0, 0.0)
            sumterm = math.complex(1.0, 0.0)
            previousterm = term
            while (math.abs(term) > 1.0e-6 and math.abs(previousterm) > 1.0e-6)
                previousterm = term
                dm += 1.0
                Cm = math.complex(dm - 1.0, 0.0)
                term = math.multiply(previousterm, math.multiply(math.divide(math.add(a, Cm), math.add(b, Cm)), math.divide(z, dm)))
                sumterm = math.add(term, sumterm)
            sumterm
        #
        # Lindhard-Sørensen correction.
        #
        lindhard: (zz, aa, bb) ->
            compton = 3.05573356675e-3  # 1.18 fm / Compton wavelength
            a3 = math.exp(math.log(aa)/3.0)
            eta = zz * @ALPHA / bb
            gg = 1.0/math.sqrt(1.0 - bb*bb)
            rho = a3*compton
            prh = bb*gg*rho
            if gg < 10.0/rho or not @switches.NuclearSize
                dk = [0.0, 0.0, 0.0]
                dmk = 0
                dkm1 = 0
                n = 1
                term1 = 0
                term2 = 0
                term3 = 1
                sumterm = 0.0
                while n < 100
                    max = if n == 1 then 2 else 1
                    for i in [0..max]
                        if i == 0
                            k = n
                        if i == 1
                            k = -n-1.0
                        if i == 2
                            k = -n
                        signk = k/math.abs(k)
                        sk = math.sqrt(k*k - @ALPHA*@ALPHA*zz*zz)
                        l = if k > 0 then k else -k - 1.0
                        Cske = math.complex(sk+1.0, eta)
                        Cketag = math.complex(k, -eta/gg)
                        Cskmeta = math.complex(sk, -eta)
                        Cexir = math.sqrt(math.divide(Cketag, Cskmeta))
                        Clg = @lngamma(Cske)
                        Cpiske = math.complex(0.0, (math.pi/2.0)*(l-sk) - Clg.im)
                        Cedr = math.multiply(Cexir, math.exp(Cpiske))
                        H = 0.0
                        Ceds = math.complex(0.0, 0.0)
                        if @switches.NuclearSize
                            Cmske = math.complex(-sk+1.0, eta)
                            Cmskmeta = math.complex(-sk, -eta)
                            Cexis = math.sqrt(math.divide(Cketag, Cmskmeta))
                            Cpimske = math.complex(0.0, (math.pi/2.0)*(l+sk) - (@lngamma(Cmske)).im)
                            Ceds = math.multiply(Cexis, math.exp(Cpimske))
                            # Caar = Cske
                            # Caas = Cmske
                            Cbbr = math.complex(2.0*sk + 1.0, 0.0)
                            Cbbs = math.complex(-2.0*sk + 1.0, 0.0)
                            Czzr = math.complex(0.0, 2.0*prh)
                            Cmprh = math.complex(0.0, -prh)
                            Clamr = math.multiply(Cexir, math.multiply(math.exp(Cmprh), @hyperg(Cske, Cbbr, Czzr)))
                            Clams = math.multiply(Cexis, math.multiply(math.exp(Cmprh), @hyperg(Cmske, Cbbs, Czzr)))
                            grgs = Clamr.im / Clams.im
                            Cgrgs = @lngamma(Cbbs)
                            grgs *= math.exp( (@lngamma(Cske)).re - (@lngamma(Cmske)).re - (@lngamma(Cbbr)).re + Cgrgs.re + 2.0*sk*math.log(2.0*prh) )
                            grgs *= -1.0 if math.cos(Cgrgs.im) < 1.0
                            if math.abs(grgs) > 1.0e-9
                                frgr = math.sqrt((gg - 1.0) / (gg + 1.0))*Clamr.re/Clamr.im
                                fsgs = math.sqrt((gg - 1.0) / (gg + 1.0))*Clams.re/Clams.im
                                gz = -1.0 * signk * (rho*gg + 1.5*@ALPHA*zz)
                                z1 = -1.0 * signk * zz
                                b0 = 1.0
                                a0 = (1.0 + 2.0*math.abs(k))*b0/(rho - gz)
                                a1 = 0.5 * (gz + rho) * b0
                                an = a1
                                anm1 = a0
                                bnm1 = b0
                                asum = a0
                                bsum = b0
                                nn = 1.0
                                while (math.abs(anm1/asum) > 1.0e-6 and math.abs(bnm1/bsum) > 1.0e-6)
                                    bn = ((rho - gz)*an + @ALPHA*z1*anm1/2.0)/(2.0*nn + 2.0*math.abs(k) + 1.0)
                                    anp1 = ((gz + rho)*bn - @ALPHA*z1*bnm1/2.0)/(2.0*nn + 2.0)
                                    asum += an
                                    bsum += bn
                                    nn += 1.0
                                    anm1 = an
                                    an = anp1
                                    bnm1 = bn
                                figi = if  k > 0 then asum/bsum else bsum/asum
                                H = ((frgr - figi)/(figi - fsgs))*grgs;
                            else
                                H = 0.0
                        #
                        # End NuclearSize
                        #
                        dk[i] = math.arg(math.add(Cedr, math.multiply(Ceds, H)))
                    if n > 1
                        dk[2] = dmk
                    sdm2 = math.sin(dk[2] - dk[1])
                    term1 = n*(n+1.0)*sdm2*sdm2/(eta*eta*(2.0*n + 1.0))
                    if n > 1
                        sd2 = math.sin(dk[0] - dkm1)
                        term1 += n*(n-1.0)*sd2*sd2/(eta*eta*(2.0*n - 1.0))
                    sdd = math.sin(dk[0] - dk[2])
                    term2 = n*sdd*sdd/(eta*eta*(4.0*n*n - 1.0))
                    term3 = term1 - 1.0/k
                    sumterm += term2 + term3
                    n += 1
                    dkm1 = dk[0]
                    dmk = dk[1]
            else
                sumterm = -math.log(prh) - 0.2  # Asymptotic value of the LS correction.
            sumterm + 0.5*bb*bb
        #
        # dE/dx
        #
        dedx: (e1, z0, a1, target) ->
            t = @absorber target
            g = 1.0 + e1/@ATOMICMASSUNIT
            delt = if @switches.NewDelta then @delta(g, t) else @olddelta(g, t)
            b2 = 1.0 - 1.0/(g*g)
            b = math.sqrt(b2)
            z1 = @effective_charge(z0, e1, t.z2)
            f1 = 0.3070722*z1*z1*t.z2/(b2*a1*t.a2)
            f2 = math.log(2.0*@ELECTRONMASS*b2/t.iadj)
            if @switches.Shell
                etam2 = 1.0/(b*b*g*g)
                cadj = 1.0e-6*(t.iadj)*(t.iadj)*etam2*(0.422377+etam2*(0.0304043-etam2*0.00038106))+1.0e-9*(t.iadj)*(t.iadj)*(t.iadj)*etam2*(3.858019+etam2*(-0.1667989 + etam2*0.00157955))
                f2 -= cadj/t.z2
                if @switches.Leung
                    f2 -= (5.0/3.0)*math.log(2.0*@ELECTRONMASS*b2/t.iadj)*(1.0e+03*t.bind/(t.z2*@ELECTRONMASS))-(t.iadj*t.iadj/(4.0*@ELECTRONMASS*@ELECTRONMASS*b2))
            f6 = 2.0*math.log(g) - b2
            #
            # Lindhard
            #
            f3 = @lindhard(z1, a1, b)
            #
            # Barkas
            #
            f4 = 1.0
            if @switches.Barkas
                v = b*g/(@ALPHA*math.sqrt(t.z2))
                fv = @fva[9]*math.exp(-2.0*math.log(v/10.0))
                f4 = 1.0 + 2.0*z1*fv/(math.sqrt(t.z2))
            # Kinematic correction
            #
            f8 = if @switches.Kinematic then 0.5*(-math.log(1.0+2.0*((5.4858e-04)*g/a1)) - ((5.4858e-04)*g/a1)*b2/(g*g)) else 0.0
            #
            # Radiative correction
            #
            f9 = if @switches.Radiative then (@ALPHA/math.pi)*b2*(6.0822 + math.log(2.0*g)*(math.log(2.0*g)*(2.4167 + 0.3333*math.log(2.0*g))-8.0314)) else 0.0
            # Placeholder for Bremsstrahlung
            Sbr = 0.0
            # Placeholder for pair production
            Spa = 0.0
            f1*(f2*f4 + f3 + f6 + (delt/2.0) + f8 + f9) + Sbr + Spa
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
    # Calculate
    #
    calculate = (eventObject) ->
        task = $('input[name=task]:checked').val()
        re = Number $('#RE').val()
        z0 = Number $('#Z').val()
        a1 = Number $('#A').val()
        target = $('#select_target').val()
        bitmask = calculate_bitmask(eventObject)
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
                type = 'd<var>E</var>/d<var>x</var>'
                result = DeDx.dedx(re, z0, a1, target)
                unit = 'A&nbsp;MeV&nbsp;g<sup>-1</sup>&nbsp;cm<sup>2</sup>'
        DeDx.nCalc += 1
        r = $('#result')
        rr = $("<tr id=\"rtr#{DeDx.nCalc}\"/>")
        $('<td/>').html(DeDx.nCalc).appendTo rr
        $('<td/>').html(type).appendTo rr
        $('<td/>').html(re).appendTo rr
        $('<td/>').html(z0).appendTo rr
        $('<td/>').html(a1).appendTo rr
        $('<td/>').html(target).appendTo rr
        $('<td/>').html(bitmask).appendTo rr
        $('<td/>').html(result).appendTo rr
        $('<td/>').html(unit).appendTo rr
        rr.appendTo r
        true
    #
    # Update bitmask
    #
    calculate_bitmask = (eventObject) ->
        bitmask = 0
        for own n, v of DeDx.switches
            vv = $("##{n}").is(':checked')
            bitmask += if vv then Number $("##{n}").val() else 0
            # console.log("#{n} -> #{vv} (was #{v})")
            DeDx.switches[n] = vv
        $('#bitmask').html(bitmask)
        bitmask
    #
    # Convert results to CSV.
    #
    csv = (eventObject) ->
        rows = []
        header = ['ID', 'Task', 'E/R', 'Z', 'A', 'Target', 'Effects', 'Result', 'Units']
        rows.push(header.join(','))
        foo = $('#result').children()
        if foo.length == 0
            alert 'No rows!'
            return false
        for row in foo
            r = []
            for col in row.children
                r.push col.innerHTML.replace(/&nbsp;/g, ' ')
            rows.push(r.join(','))
        c = rows.join('\r\n') + '\r\n'
        download = $("<a/>", {href: 'data:text/csv;charset=utf-8,' + encodeURIComponent(c), download: "crange.csv"})
        download.appendTo $("#downloadCSV")
        download[0].click()
        $("#downloadCSV").empty()
        true
    #
    # Reset the form and clear results.
    #
    resetForm = (eventObject) ->
        document.getElementById('recalc').reset()
        $('#result').empty()
        DeDx.nCalc = 0
        for own n, v of DeDx.switches
            DeDx.switches[n] = false
        DeDx.switches.NewDelta = true
        DeDx.switches.FiniteNuclearSize = true
        $('#bitmask').html(40)
        true
    #
    # Load the target data.
    #
    if DeDx.targets.length == 0
        $.getJSON("target.json",{},(data) -> DeDx.targets = data).fail(
            () -> alert("JSON error!")
            ).done(absorberTable)
    #
    # Bind to various changes.
    #
    for div in DeDx.validateDivs
        if $("##{div.name}").length > 0
            $("##{div.name}").change(div, validateNumber).change()
    $('#div_switch input[type=checkbox]').click(calculate_bitmask)
    $('#calculate').click(calculate)
    $('#resetForm').click(resetForm)
    $('#CSV').click(csv)
    true
