$(function () {
    let DeDx = {
        nCalc: 0,
        ALPHA: 7.29735301383e-3,
        ATOMICMASSUNIT: 931.4943,
        PROTONMASS: 938.2723,
        ELECTRONMASS: 0.511003e+6,
        fva: [0.33, 0.078, 0.03, 0.014, 0.0084, 0.0053, 0.0035, 0.0025, 0.0019, 0.0014],
        targets: [],
        switches: {
            Barkas: false,
            Shell: false,
            Leung: false,
            NewDelta: true,
            NewElectronCapture: false,
            FiniteNuclearSize: true,
            Kinematic: false,
            Radiative: false,
            Pair: false,
            Bremsstrahlung: false
        },
        validateDivs: [{
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
            }],
        absorber: function (target) {
            let t;
            for (let j = 0; j < this.targets.length; j++) {
                t = this.targets[j];
                if (t.name === target)
                    return t;
            }
            return t;
        },
        effective_charge: function (z0, e1, z2) {
            const g = 1.0 + e1 / this.ATOMICMASSUNIT;
            const b2 = 1.0 - 1.0 / (g * g);
            const b = math.sqrt(b2);
            let f1;
            if (this.switches.NewElectronCapture) {
                if (z2 === 4.0) {
                    f1 = (2.045 + 2.000 * math.exp(-0.04369 * z0)) * math.exp(-7.000 * math.exp(0.2643 * math.log(e1)) * math.exp(-0.4171 * math.log(z0)));
                }
                else {
                    if (z2 === 6.0) {
                        f1 = (2.584 + 1.910 * math.exp(-0.03958 * z0)) * math.exp(-6.933 * math.exp(0.2433 * math.log(e1)) * math.exp(-0.3969 * math.log(z0)));
                    }
                    else {
                        f1 = ((1.164 + 0.2319 * math.exp(-0.004302 * z2)) + 1.658 * math.exp(-0.05170 * z0)) * math.exp(-(8.144 + 0.9876 * math.log(z2)) * math.exp((0.3140 + 0.01072 * math.log(z2)) * math.log(e1)) * math.exp(-(0.5218 + 0.02521 * math.log(z2)) * math.log(z0)));
                    }
                }
            }
            else {
                let z23 = math.exp((2.0 / 3.0) * math.log(z0));
                let capA = 1.16 - z2 * (1.91e-03 - 1.26e-05 * z2);
                let capB = (1.18 - z2 * (7.5e-03 - 4.53e-05 * z2)) / this.ALPHA;
                f1 = capA * math.exp(-capB * b / z23);
            }
            return (1.0 - f1) * z0;
        },
        delta: function (g, t) {
            let X0 = t.X0;
            let X1 = t.X1;
            let cbar = 2.0 * math.log(t.iadj / t.pla) + 1.0;
            const b = math.sqrt(1.0 - 1.0 / (g * g));
            const X = math.log10(b * g);
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
        olddelta: function (g, t) {
            if (g < 1.8)
                return 0.0;
            const cbar = 2.0 * math.log(t.iadj / t.pla) + 1.0;
            const b = math.sqrt(1.0 - 1.0 / (g * g));
            let y = 2.0 * math.log(b * g);
            let y0;
            let y1;
            if (t.etad > 0) {
                y += math.log(t.etad);
                if (cbar >= 12.25) {
                    y1 = 23.03;
                    y0 = cbar >= 13.804 ? 1.502 * cbar - 11.52 : 9.212;
                }
                else {
                    y1 = 18.42;
                    if (cbar < 12.25)
                        y0 = 9.212;
                    if (cbar < 11.5)
                        y0 = 8.751;
                    if (cbar < 11.0)
                        y0 = 8.291;
                    if (cbar < 10.5)
                        y0 = 7.830;
                    if (cbar < 10.0)
                        y0 = 7.370;
                }
            }
            else {
                if (t.iadj >= 100.0) {
                    y1 = 13.82;
                    y0 = cbar >= 5.215 ? 1.502 * cbar - 6.909 : 0.9212;
                }
                else {
                    y1 = 9.212;
                    y0 = cbar >= 3.681 ? 1.502 * cbar - 4.606 : 0.9212;
                }
            }
            if (y < y0) {
                return 0.0;
            }
            else {
                if (y > y1) {
                    return y - cbar;
                }
                else {
                    let dy3 = (y1 - y0) * (y1 - y0) * (y1 - y0);
                    let a = (cbar - y0) / dy3;
                    return y - cbar + a * (y1 - y) * (y1 - y) * (y1 - y);
                }
            }
        },
        lngamma: function (z) {
            const coeff = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
            let x;
            let y;
            if (z.re > 0) {
                x = z.re - 1.0;
                y = z.im;
            }
            else {
                x = -z.re;
                y = -z.im;
            }
            let r = math.sqrt((x + 5.5) * (x + 5.5) + y * y);
            let aterm1 = y * math.log(r);
            let aterm2 = (x + 0.5) * math.atan2(y, x + 5.5) - y;
            let lterm1 = (x + 0.5) * math.log(r);
            let lterm2 = -y * math.atan2(y, x + 5.5) - (x + 5.5) + 0.5 * math.log(2.0 * math.pi);
            let num = 0.0;
            let denom = 1.000000000190015;
            for (let i = 0; i < coeff.length; i++) {
                let cterm = coeff[i] / ((x + i + 1) * (x + i + 1) + y * y);
                num += cterm;
                denom += (x + i + 1) * cterm;
            }
            num *= -y;
            let aterm3 = math.atan2(num, denom);
            let lterm3 = 0.5 * math.log(num * num + denom * denom);
            let result = math.complex(lterm1 + lterm2 + lterm3, aterm1 + aterm2 + aterm3);
            if (z.re < 0) {
                let lpi = math.complex(math.log(math.pi), 0);
                result = math.subtract(lpi, math.add(result, math.log(math.sin(math.multiply(math.pi, z)))));
            }
            return result;
        },
        hyperg: function (a, b, z) {
            let dm = 0.0;
            let term = math.complex(1.0, 0.0);
            let sumterm = math.complex(1.0, 0.0);
            let previousterm = term;
            while (math.abs(term) > 1.0e-6 && math.abs(previousterm) > 1.0e-6) {
                previousterm = term;
                dm += 1.0;
                let Cm = math.complex(dm - 1.0, 0.0);
                term = math.multiply(previousterm, math.multiply(math.divide(math.add(a, Cm), math.add(b, Cm)), math.divide(z, dm)));
                sumterm = math.add(term, sumterm);
            }
            return sumterm;
        },
        lindhard: function (zz, aa, bb) {
            const compton = 3.05573356675e-3;
            const a3 = math.exp(math.log(aa) / 3.0);
            const eta = zz * this.ALPHA / bb;
            const gg = 1.0 / math.sqrt(1.0 - bb * bb);
            const rho = a3 * compton;
            const prh = bb * gg * rho;
            let n = 1;
            let sumterm = 0.0;
            let term1 = 0;
            let term3 = 1;
            let term2 = 0;
            if (gg < 10.0 / rho || !this.switches.FiniteNuclearSize) {
                let dk = [0.0, 0.0, 0.0];
                let dmk = 0;
                let dkm1 = 0;
                while (n < 100) {
                    let max = n === 1 ? 3 : 2;
                    for (let i = 0; i < max; i++) {
                        let k;
                        if (i === 0)
                            k = n;
                        if (i === 1)
                            k = -n - 1.0;
                        if (i === 2)
                            k = -n;
                        let signk = k / math.abs(k);
                        let sk = math.sqrt(k * k - this.ALPHA * this.ALPHA * zz * zz);
                        let l = k > 0 ? k : -k - 1.0;
                        let Cske = math.complex(sk + 1.0, eta);
                        let Cketag = math.complex(k, -eta / gg);
                        let Cskmeta = math.complex(sk, -eta);
                        let Cexir = math.sqrt(math.divide(Cketag, Cskmeta));
                        let Clg = this.lngamma(Cske);
                        let Cpiske = math.complex(0.0, (math.pi / 2.0) * (l - sk) - Clg.im);
                        let Cedr = math.multiply(Cexir, math.exp(Cpiske));
                        let H = 0.0;
                        let Ceds = math.complex(0.0, 0.0);
                        if (this.switches.FiniteNuclearSize) {
                            let Cmske = math.complex(-sk + 1.0, eta);
                            let Cmskmeta = math.complex(-sk, -eta);
                            let Cexis = math.sqrt(math.divide(Cketag, Cmskmeta));
                            let Cpimske = math.complex(0.0, (math.pi / 2.0) * (l + sk) - (this.lngamma(Cmske)).im);
                            Ceds = math.multiply(Cexis, math.exp(Cpimske));
                            let Cbbr = math.complex(2.0 * sk + 1.0, 0.0);
                            let Cbbs = math.complex(-2.0 * sk + 1.0, 0.0);
                            let Czzr = math.complex(0.0, 2.0 * prh);
                            let Cmprh = math.complex(0.0, -prh);
                            let Clamr = math.multiply(Cexir, math.multiply(math.exp(Cmprh), this.hyperg(Cske, Cbbr, Czzr)));
                            let Clams = math.multiply(Cexis, math.multiply(math.exp(Cmprh), this.hyperg(Cmske, Cbbs, Czzr)));
                            let grgs = Clamr.im / Clams.im;
                            let Cgrgs = this.lngamma(Cbbs);
                            grgs *= math.exp((this.lngamma(Cske)).re - (this.lngamma(Cmske)).re - (this.lngamma(Cbbr)).re + Cgrgs.re + 2.0 * sk * math.log(2.0 * prh));
                            if (math.cos(Cgrgs.im) < 1.0)
                                grgs *= -1.0;
                            if (math.abs(grgs) > 1.0e-9) {
                                let frgr = math.sqrt((gg - 1.0) / (gg + 1.0)) * Clamr.re / Clamr.im;
                                let fsgs = math.sqrt((gg - 1.0) / (gg + 1.0)) * Clams.re / Clams.im;
                                let gz = -1.0 * signk * (rho * gg + 1.5 * this.ALPHA * zz);
                                let z1 = -1.0 * signk * zz;
                                let b0 = 1.0;
                                let a0 = (1.0 + 2.0 * math.abs(k)) * b0 / (rho - gz);
                                let a1 = 0.5 * (gz + rho) * b0;
                                let an = a1;
                                let anm1 = a0;
                                let bnm1 = b0;
                                let asum = a0;
                                let bsum = b0;
                                let nn = 1.0;
                                while (math.abs(anm1 / asum) > 1.0e-6 && math.abs(bnm1 / bsum) > 1.0e-6) {
                                    let bn = ((rho - gz) * an + this.ALPHA * z1 * anm1 / 2.0) / (2.0 * nn + 2.0 * math.abs(k) + 1.0);
                                    let anp1 = ((gz + rho) * bn - this.ALPHA * z1 * bnm1 / 2.0) / (2.0 * nn + 2.0);
                                    asum += an;
                                    bsum += bn;
                                    nn += 1.0;
                                    anm1 = an;
                                    an = anp1;
                                    bnm1 = bn;
                                }
                                let figi = k > 0 ? asum / bsum : bsum / asum;
                                H = ((frgr - figi) / (figi - fsgs)) * grgs;
                            }
                            else {
                                H = 0.0;
                            }
                        }
                        dk[i] = math.arg(math.add(Cedr, math.multiply(Ceds, H)));
                    }
                    if (n > 1)
                        dk[2] = dmk;
                    let sdm2 = math.sin(dk[2] - dk[1]);
                    term1 = n * (n + 1.0) * sdm2 * sdm2 / (eta * eta * (2.0 * n + 1.0));
                    if (n > 1) {
                        let sd2 = math.sin(dk[0] - dkm1);
                        term1 += n * (n - 1.0) * sd2 * sd2 / (eta * eta * (2.0 * n - 1.0));
                    }
                    let sdd = math.sin(dk[0] - dk[2]);
                    term2 = n * sdd * sdd / (eta * eta * (4.0 * n * n - 1.0));
                    term3 = term1 - 1.0 / n;
                    sumterm += term2 + term3;
                    n += 1;
                    dkm1 = dk[0];
                    dmk = dk[1];
                }
            }
            else {
                sumterm = -math.log(prh) - 0.2;
            }
            return sumterm + 0.5 * bb * bb;
        },
        stop: function (e1, z0, a1, target) {
            let t = this.absorber(target);
            const g = 1.0 + e1 / this.ATOMICMASSUNIT;
            const delt = this.switches.NewDelta ? this.delta(g, t) : this.olddelta(g, t);
            const b2 = 1.0 - 1.0 / (g * g);
            const b = math.sqrt(b2);
            const z1 = this.effective_charge(z0, e1, t.z2);
            let f1 = 0.3070722 * z1 * z1 * t.z2 / (b2 * a1 * t.a2);
            let f2 = math.log(2.0 * this.ELECTRONMASS * b2 / t.iadj);
            if (this.switches.Shell) {
                let etam2 = 1.0 / (b * b * g * g);
                let cadj = 1.0e-6 * t.iadj * t.iadj * etam2 * (0.422377 + etam2 * (0.0304043 - etam2 * 0.00038106)) + 1.0e-9 * t.iadj * t.iadj * t.iadj * etam2 * (3.858019 + etam2 * (-0.1667989 + etam2 * 0.00157955));
                f2 -= cadj / t.z2;
                if (this.switches.Leung) {
                    f2 -= (5.0 / 3.0) * math.log(2.0 * this.ELECTRONMASS * b2 / t.iadj) * (1.0e+03 * t.bind / (t.z2 * this.ELECTRONMASS)) - (t.iadj * t.iadj / (4.0 * this.ELECTRONMASS * this.ELECTRONMASS * b2));
                }
            }
            let f6 = 2.0 * math.log(g) - b2;
            let f3 = this.lindhard(z1, a1, b);
            let f4 = 1.0;
            if (this.switches.Barkas) {
                let v = b * g / (this.ALPHA * math.sqrt(t.z2));
                let fv = this.fva[9] * math.exp(-2.0 * math.log(v / 10.0));
                f4 = 1.0 + 2.0 * z1 * fv / (math.sqrt(t.z2));
            }
            let f8 = this.switches.Kinematic ? 0.5 * (-math.log(1.0 + 2.0 * (5.4858e-04 * g / a1)) - (5.4858e-04 * g / a1) * b2 / (g * g)) : 0.0;
            let f9 = this.switches.Radiative ? (this.ALPHA / math.pi) * b2 * (6.0822 + math.log(2.0 * g) * (math.log(2.0 * g) * (2.4167 + 0.3333 * math.log(2.0 * g)) - 8.0314)) : 0.0;
            let Spa = 0.0;
            if (this.switches.Pair) {
                let dpa = 1.0 / math.sqrt(g);
                let ldpa = math.log(dpa);
                let l0 = math.log(2.0 * g);
                let Lpa0s = (19.0 / 9.0) * math.log(183.0 * math.exp(-1.0 / 3.0 * math.log(t.z2)) / (1.0 + 4.0 * 6.25470095193633 * 183.0 * math.exp(-1.0 / 3.0 * math.log(t.z2)) / g));
                let Lpa1 = dpa * (4178.0 / (81 * math.pi * math.pi) - 21.0 / 27.0 - 248.0 * l0 / (27.0 * math.pi * math.pi)
                    + (28.0 * l0 / 9.0 - 446.0 / 27.0) * 2.0 * ldpa / (math.pi * math.pi) + 14.0 * 4.0 * ldpa * ldpa / (9.0 * math.pi * math.pi));
                let Lpa = Lpa0s + Lpa1;
                Spa = 4.08803936906434e-06 * (z1 * z1 / a1) * (t.z2 * t.z2 / t.a2) * (1.0 + 1.0 / t.z2) * g * Lpa;
            }
            let Sbr = 0.0;
            if (this.switches.Bremsstrahlung) {
                let Bbr = math.log(1.0 + 2.0 * g * 0.179524783764566 / (math.exp((1.0 / 3.0) * math.log(a1)) + math.exp((1.0 / 3.0) * math.log(t.a2))) / a1);
                Sbr = 5.21721169334564e-07 * (z1 * z1 / a1) * (z1 * z1 / a1) * (t.z2 * t.z2 / t.a2) * g * Bbr;
            }
            return f1 * (f2 * f4 + f3 + f6 - (delt / 2.0) + f8 + f9) + Sbr + Spa;
        }
    };
    let validateNumber = function (eventObject) {
        if ($('#previous_target').val() !== 'Unknown')
            eventObject.data.first = false;
        let div = eventObject.data.name;
        if (eventObject.data.type === 'radio') {
            let name = div.split('_')[1];
            eventObject.data.valid = $("input[name=" + name + "]:checked").length === 1;
        }
        else {
            let input = $('#' + div.split('_')[1].toUpperCase());
            let patt = eventObject.data.type === 'int' ? new RegExp(/^[0-9]+$/i) : new RegExp(/^[0-9]+(\.[0-9]*|)(e[+-]?[0-9]+|)$/i);
            eventObject.data.valid = patt.test(input.val().toString());
            if (!eventObject.data.first) {
                if (eventObject.data.valid) {
                    input.removeClass('is-invalid').addClass('is-valid');
                    $('#' + div + '_helpblock').html(eventObject.data.help);
                }
                else {
                    input.removeClass('is-valid').addClass('is-invalid');
                    $('#' + div + '_helpblock').html(eventObject.data.help + " " + (eventObject.data.type === 'int' ? 'Integer' : 'Float') + " value required!");
                }
            }
        }
        eventObject.data.first = false;
        let validity = [];
        for (let j = 0; j < DeDx.validateDivs.length; j++) {
            validity.push(DeDx.validateDivs[j].valid);
        }
        let formvalid = validity.every(function (currentValue) { return currentValue; });
        $('#calculate').prop('disabled', !formvalid);
        return eventObject.data.valid;
    };
    let absorberTable = function () {
        let select = $('#select_target');
        let previous_target = $('#previous_target');
        let selected_target = previous_target.length === 1 ? previous_target.val().toString() : 'Unknown';
        let k = 0;
        for (let j = 0; j < DeDx.targets.length; j++) {
            if (DeDx.targets[j].name !== 'Unknown') {
                let rowid = 't' + k;
                k += 1;
                $("<option id=\"" + rowid + "\"/>").html(DeDx.targets[j].name).appendTo(select);
                if (DeDx.targets[j].name === selected_target) {
                    $("#" + rowid).prop('selected', true);
                }
            }
        }
    };
    let calculate = function (eventObject) {
        let task = $('input[name=task]:checked').val().toString();
        let re = Number($('#RE').val());
        let z0 = Number($('#Z').val());
        let a1 = Number($('#A').val());
        let target = $('#select_target').val().toString();
        let bitmask = calculate_bitmask(eventObject);
        let type;
        let result;
        let unit;
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
                result = DeDx.stop(re, z0, a1, target);
                unit = 'A&nbsp;MeV&nbsp;g<sup>-1</sup>&nbsp;cm<sup>2</sup>';
        }
        DeDx.nCalc += 1;
        let r = $('#result');
        let rr = $("<tr id=\"rtr" + DeDx.nCalc + "\"/>");
        $('<td/>').html(DeDx.nCalc.toString()).appendTo(rr);
        $('<td/>').html(type).appendTo(rr);
        $('<td/>').html(re.toString()).appendTo(rr);
        $('<td/>').html(z0.toString()).appendTo(rr);
        $('<td/>').html(a1.toString()).appendTo(rr);
        $('<td/>').html(target).appendTo(rr);
        $('<td/>').html(bitmask.toString()).appendTo(rr);
        $('<td/>').html(result.toString()).appendTo(rr);
        $('<td/>').html(unit).appendTo(rr);
        rr.appendTo(r);
    };
    let calculate_bitmask = function (_eventObject) {
        let bitmask = 0;
        let n;
        for (n in DeDx.switches) {
            if (!DeDx.switches.hasOwnProperty(n))
                continue;
            let vv = $("#" + n).is(':checked');
            bitmask += vv ? Number($("#" + n).val()) : 0;
            DeDx.switches[n] = vv;
        }
        $('#bitmask').html(bitmask.toString());
        return bitmask;
    };
    let resetForm = function (_eventObject) {
        let recalc = document.getElementById('recalc');
        recalc.reset();
        $('#result').empty();
        DeDx.nCalc = 0;
        for (let j = 0; j < DeDx.validateDivs.length; j++) {
            if (DeDx.validateDivs[j].type !== "radio") {
                let name = DeDx.validateDivs[j].name.split('_')[1].toUpperCase();
                $("#" + name).removeClass('is-invalid').removeClass('is-valid');
                $("#" + DeDx.validateDivs[j].name + "_helpblock").html(DeDx.validateDivs[j].help);
            }
        }
        let n;
        for (n in DeDx.switches) {
            if (!DeDx.switches.hasOwnProperty(n))
                continue;
            DeDx.switches[n] = false;
        }
        DeDx.switches.NewDelta = true;
        DeDx.switches.FiniteNuclearSize = true;
        $('#bitmask').html("40");
    };
    let csv = function (_eventObject) {
        let rows = [];
        const header = ['ID', 'Task', 'E/R', 'Z', 'A', 'Target', 'Effects', 'Result', 'Units'];
        rows.push(header.join(','));
        let foo = $('#result').children();
        if (foo.length === 0) {
            alert('No rows!');
            return;
        }
        for (let j = 0; j < foo.length; j++) {
            let r = [];
            for (let m = 0; m < foo[j].children.length; m++) {
                r.push(foo[j].children[m].innerHTML.replace(/&nbsp;/g, ' '));
            }
            rows.push(r.join(','));
        }
        let c = rows.join('\r\n') + '\r\n';
        let download = $("<a/>", {
            href: 'data:text/csv;charset=utf-8,' + encodeURIComponent(c),
            download: "crange.csv"
        });
        download.appendTo($("#downloadCSV"));
        download[0].click();
        $("#downloadCSV").empty();
    };
    if (DeDx.targets.length === 0) {
        $.getJSON("target.json", {}, function (data) {
            return DeDx.targets = data;
        }).fail(function () {
            return alert("JSON error!");
        }).done(absorberTable);
    }
    for (let j = 0; j < DeDx.validateDivs.length; j++) {
        let div = DeDx.validateDivs[j];
        if ($("#" + div.name).length > 0) {
            $("#" + div.name).change(div, validateNumber).change();
        }
    }
    $('input[type=checkbox]').click(calculate_bitmask);
    $('#calculate').click(calculate);
    $('#resetForm').click(resetForm);
    $('#CSV').click(csv);
});
