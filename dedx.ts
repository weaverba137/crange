/// <reference path="../../../../local/products/node_modules/lib/node_modules/@types/jquery/index.d.ts" />
/// <reference path="../../../../local/products/node_modules/lib/node_modules/@types/mathjs/index.d.ts" />
type FVA = [number, number, number, number, number, number, number, number, number, number];

interface Switch {
    Barkas: boolean;
    Shell: boolean;
    Leung: boolean;
    NewDelta: boolean;
    NewElectronCapture: boolean;
    FiniteNuclearSize: boolean;
    Kinematic: boolean;
    Radiative: boolean;
    Pair: boolean;
    Bremsstrahlung: boolean;
}

interface Validation {
    name: string;
    type: string;
    valid: boolean;
    first: boolean;
    help: string;
}

interface DEDX {
    nCalc: number;
    ALPHA: number;
    ATOMICMASSUNIT: number;
    PROTONMASS: number;
    ELECTRONMASS: number;
    fva: FVA;
    targets: Target[];
    switches: Switch;
    validateDivs: Validation[];
    absorber(target: string): Target;
    effective_charge(z0: number, e1: number, z2: number): number;
    delta(g: number, t: Target): number;
    olddelta(g: number, t: Target): number;
    lngamma(z: math.Complex): math.Complex;
    hyperg(a: math.Complex, b: math.Complex, z: math.Complex): math.Complex;
    lindhard(zz: number, aa: number, bb: number): number;
    stop(e1: number, z0: number, a1: number, target: string): number;
}

declare namespace JQuery {
    interface Event {
        data: Validation;
    }
}

$(
    function(): void {
        let DeDx: DEDX = {
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
            absorber: function(target: string): Target {
                let t: Target;
                for (let j: number = 0; j < this.targets.length; j++) {
                    t = this.targets[j];
                    // console.log("t.name = '" + t.name + "'");
                    if (t.name === target) return t;
                }
                // alert("Target "+target+" not found!");
                return t;
            },
            effective_charge: function(z0: number, e1: number, z2: number): number {
                const g: number = 1.0 + e1 / this.ATOMICMASSUNIT;
                const b2: number = 1.0 - 1.0 / (g * g);
                const b: number = math.sqrt(b2);
                let f1: number;
                if (this.switches.NewElectronCapture) {
                    if (z2 === 4.0) {
                        f1 = (2.045 + 2.000 * math.exp(-0.04369 * z0)) * math.exp(-7.000 * math.exp(0.2643 * math.log(e1)) * math.exp(-0.4171 * math.log(z0)));
                    } else {
                        if (z2 === 6.0) {
                            f1 = (2.584 + 1.910 * math.exp(-0.03958 * z0)) * math.exp(-6.933 * math.exp(0.2433 * math.log(e1)) * math.exp(-0.3969 * math.log(z0)));
                        } else {
                            f1 = ((1.164 + 0.2319 * math.exp(-0.004302 * z2)) + 1.658 * math.exp(-0.05170 * z0)) * math.exp(-(8.144 + 0.9876 * math.log(z2)) * math.exp((0.3140 + 0.01072 * math.log(z2)) * math.log(e1)) * math.exp(-(0.5218 + 0.02521 * math.log(z2)) * math.log(z0)));
                        }
                    }
                } else {
                    let z23: number = math.exp((2.0 / 3.0) * math.log(z0));
                    // const capA: number = 1.0;
                    let capA: number = 1.16 - z2 * (1.91e-03 - 1.26e-05 * z2);
                    // const capB: number = 130.0;
                    let capB: number = (1.18 - z2 * (7.5e-03 - 4.53e-05 * z2)) / this.ALPHA;
                    f1 = capA * math.exp(-capB * b / z23);
                }
                return (1.0 - f1) * z0;
            },
            delta: function(g: number, t: Target): number {
                let X0: number = t.X0;
                let X1: number = t.X1;
                let cbar: number = 2.0 * math.log(t.iadj / t.pla) + 1.0;
                const b: number = math.sqrt(1.0 - 1.0 / (g * g));
                const X: number = math.log10(b * g);
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
            olddelta: function(g: number, t: Target): number {
                if (g < 1.8) return 0.0;
                const cbar: number = 2.0 * math.log(t.iadj / t.pla) + 1.0;
                const b: number = math.sqrt(1.0 - 1.0 / (g * g));
                let y: number = 2.0 * math.log(b * g);
                let y0: number;
                let y1: number;
                if (t.etad > 0) {
                    y += math.log(t.etad);
                    if (cbar >= 12.25) {
                        y1 = 23.03;
                        y0 = cbar >= 13.804 ? 1.502 * cbar - 11.52 : 9.212;
                    } else {
                        y1 = 18.42;
                        if (cbar < 12.25) y0 = 9.212;
                        if (cbar < 11.5) y0 = 8.751;
                        if (cbar < 11.0) y0 = 8.291;
                        if (cbar < 10.5) y0 = 7.830;
                        if (cbar < 10.0) y0 = 7.370;
                    }
                } else {
                    if (t.iadj >= 100.0) {
                        y1 = 13.82;
                        y0 = cbar >= 5.215 ? 1.502 * cbar - 6.909 : 0.9212;
                    } else {
                        y1 = 9.212;
                        y0 = cbar >= 3.681 ? 1.502 * cbar - 4.606 : 0.9212;
                    }
                }
                if (y < y0) {
                    return 0.0;
                } else {
                    if (y > y1) {
                        return y - cbar;
                    } else {
                        let dy3: number = (y1 - y0) * (y1 - y0) * (y1 - y0);
                        let a: number = (cbar - y0) / dy3;
                        return y - cbar + a * (y1 - y) * (y1 - y) * (y1 - y);
                    }
                }
            },
            lngamma: function(z: math.Complex): math.Complex {
                const coeff: number[] = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
                let x: number;
                let y: number;
                if (z.re > 0) {
                    x = z.re - 1.0;
                    y = z.im;
                } else {
                    x = -z.re;
                    y = -z.im;
                }
                let r: number = math.sqrt((x + 5.5) * (x + 5.5) + y * y);
                let aterm1: number = y * math.log(r);
                let aterm2: number = (x + 0.5) * math.atan2(y, x + 5.5) - y;
                let lterm1: number = (x + 0.5) * math.log(r);
                let lterm2: number = -y * math.atan2(y, x + 5.5) - (x + 5.5) + 0.5 * math.log(2.0 * math.pi);
                let num: number = 0.0;
                let denom: number = 1.000000000190015;
                for (let i: number = 0; i < coeff.length; i++) {
                    let cterm: number = coeff[i] / ((x + i + 1) * (x + i + 1) + y * y);
                    num += cterm;
                    denom += (x + i + 1) * cterm;
                }
                num *= -y;
                let aterm3: number = math.atan2(num, denom);
                let lterm3: number = 0.5 * math.log(num * num + denom * denom);
                let result: math.Complex = math.complex(lterm1 + lterm2 + lterm3, aterm1 + aterm2 + aterm3);
                if (z.re < 0) {
                    let lpi: math.Complex = math.complex(math.log(math.pi), 0);
                    result = math.subtract(lpi, math.add(result, math.log(math.sin(math.multiply(math.pi, z))))) as math.Complex;
                }
                return result;
            },
            hyperg: function(a: math.Complex, b: math.Complex, z: math.Complex): math.Complex {
                let dm: number = 0.0;
                let term: math.Complex = math.complex(1.0, 0.0);
                let sumterm: math.Complex = math.complex(1.0, 0.0);
                let previousterm: math.Complex = term;
                while (math.abs(term) > 1.0e-6 && math.abs(previousterm) > 1.0e-6) {
                    previousterm = term;
                    dm += 1.0;
                    let Cm: math.Complex = math.complex(dm - 1.0, 0.0);
                    term = math.multiply(previousterm, math.multiply(math.divide(math.add(a, Cm), math.add(b, Cm)), math.divide(z, dm))) as math.Complex;
                    sumterm = math.add(term, sumterm) as math.Complex;
                }
                return sumterm;
            },
            lindhard: function(zz: number, aa: number, bb: number): number {
                const compton: number = 3.05573356675e-3;
                const a3: number = math.exp(math.log(aa) / 3.0);
                const eta: number = zz * this.ALPHA / bb;
                const gg: number = 1.0 / math.sqrt(1.0 - bb * bb);
                const rho: number = a3 * compton;
                const prh: number = bb * gg * rho;
                let n: number = 1;
                let sumterm: number = 0.0;
                let term1: number = 0;
                let term3: number = 1;
                let term2: number = 0;
                if (gg < 10.0 / rho || !this.switches.FiniteNuclearSize) {
                    let dk: number[] = [0.0, 0.0, 0.0];
                    let dmk: number = 0;
                    let dkm1: number = 0;
                    while (n < 100) {
                        let max: number = n === 1 ? 3 : 2;
                        for (let i: number = 0; i < max; i++) {
                            let k: number;
                            if (i === 0) k = n;
                            if (i === 1) k = -n - 1.0;
                            if (i === 2) k = -n;
                            let signk: number = k / math.abs(k);
                            let sk: number = math.sqrt(k * k - this.ALPHA * this.ALPHA * zz * zz);
                            let l: number = k > 0 ? k : -k - 1.0;
                            let Cske: math.Complex = math.complex(sk + 1.0, eta);
                            let Cketag: math.Complex = math.complex(k, -eta / gg);
                            let Cskmeta: math.Complex = math.complex(sk, -eta);
                            let Cexir: math.Complex = math.sqrt(math.divide(Cketag, Cskmeta)) as math.Complex;
                            let Clg: math.Complex = this.lngamma(Cske);
                            let Cpiske: math.Complex = math.complex(0.0, (math.pi / 2.0) * (l - sk) - Clg.im);
                            let Cedr: math.Complex = math.multiply(Cexir, math.exp(Cpiske)) as math.Complex;
                            let H: number = 0.0;
                            let Ceds: math.Complex = math.complex(0.0, 0.0);
                            if (this.switches.FiniteNuclearSize) {
                                let Cmske: math.Complex = math.complex(-sk + 1.0, eta);
                                let Cmskmeta: math.Complex = math.complex(-sk, -eta);
                                let Cexis: math.Complex = math.sqrt(math.divide(Cketag, Cmskmeta)) as math.Complex;
                                let Cpimske: math.Complex = math.complex(0.0, (math.pi / 2.0) * (l + sk) - (this.lngamma(Cmske)).im);
                                Ceds = math.multiply(Cexis, math.exp(Cpimske)) as math.Complex;
                                let Cbbr: math.Complex = math.complex(2.0 * sk + 1.0, 0.0);
                                let Cbbs: math.Complex = math.complex(-2.0 * sk + 1.0, 0.0);
                                let Czzr: math.Complex = math.complex(0.0, 2.0 * prh);
                                let Cmprh: math.Complex = math.complex(0.0, -prh);
                                let Clamr: math.Complex = math.multiply(Cexir, math.multiply(math.exp(Cmprh), this.hyperg(Cske, Cbbr, Czzr))) as math.Complex;
                                let Clams: math.Complex = math.multiply(Cexis, math.multiply(math.exp(Cmprh), this.hyperg(Cmske, Cbbs, Czzr))) as math.Complex;
                                let grgs: number  = Clamr.im / Clams.im;
                                let Cgrgs: math.Complex = this.lngamma(Cbbs);
                                grgs *= math.exp((this.lngamma(Cske)).re - (this.lngamma(Cmske)).re - (this.lngamma(Cbbr)).re + Cgrgs.re + 2.0 * sk * math.log(2.0 * prh));
                                if (math.cos(Cgrgs.im) < 1.0) grgs *= -1.0;
                                if (math.abs(grgs) > 1.0e-9) {
                                    let frgr: number = math.sqrt((gg - 1.0) / (gg + 1.0)) * Clamr.re / Clamr.im;
                                    let fsgs: number = math.sqrt((gg - 1.0) / (gg + 1.0)) * Clams.re / Clams.im;
                                    let gz: number = -1.0 * signk * (rho * gg + 1.5 * this.ALPHA * zz);
                                    let z1: number = -1.0 * signk * zz;
                                    let b0: number = 1.0;
                                    let a0: number = (1.0 + 2.0 * math.abs(k)) * b0 / (rho - gz);
                                    let a1: number = 0.5 * (gz + rho) * b0;
                                    let an: number = a1;
                                    let anm1: number = a0;
                                    let bnm1: number = b0;
                                    let asum: number = a0;
                                    let bsum: number = b0;
                                    let nn: number = 1.0;
                                    while (math.abs(anm1 / asum) > 1.0e-6 && math.abs(bnm1 / bsum) > 1.0e-6) {
                                        let bn: number = ((rho - gz) * an + this.ALPHA * z1 * anm1 / 2.0) / (2.0 * nn + 2.0 * math.abs(k) + 1.0);
                                        let anp1: number = ((gz + rho) * bn - this.ALPHA * z1 * bnm1 / 2.0) / (2.0 * nn + 2.0);
                                        asum += an;
                                        bsum += bn;
                                        nn += 1.0;
                                        anm1 = an;
                                        an = anp1;
                                        bnm1 = bn;
                                    }
                                    let figi: number = k > 0 ? asum / bsum : bsum / asum;
                                    H = ((frgr - figi) / (figi - fsgs)) * grgs;
                                    // console.log("H = " + H.toString())
                                } else {
                                    H = 0.0;
                                }
                            }
                            dk[i] = math.arg(math.add(Cedr, math.multiply(Ceds, H))) as number;
                        }
                        if (n > 1) dk[2] = dmk;
                        let sdm2: number = math.sin(dk[2] - dk[1]);
                        term1 = n * (n + 1.0) * sdm2 * sdm2 / (eta * eta * (2.0 * n + 1.0));
                        if (n > 1) {
                            let sd2: number = math.sin(dk[0] - dkm1);
                            term1 += n * (n - 1.0) * sd2 * sd2 / (eta * eta * (2.0 * n - 1.0));
                        }
                        let sdd: number = math.sin(dk[0] - dk[2]);
                        term2 = n * sdd * sdd / (eta * eta * (4.0 * n * n - 1.0));
                        term3 = term1 - 1.0 / n;
                        sumterm += term2 + term3;
                        n += 1;
                        dkm1 = dk[0];
                        dmk = dk[1];
                    }
                } else {
                    sumterm = -math.log(prh) - 0.2;
                }
                return sumterm + 0.5 * bb * bb;
            },
            stop: function(e1: number, z0: number, a1: number, target: string): number {
                let t: Target = this.absorber(target);
                // console.log("t.name = " + t.name);
                const g: number = 1.0 + e1 / this.ATOMICMASSUNIT;
                const delt: number = this.switches.NewDelta ? this.delta(g, t) : this.olddelta(g, t);
                // console.log("delt = " + delt.toString());
                const b2: number = 1.0 - 1.0 / (g * g);
                const b: number = math.sqrt(b2);
                const z1: number = this.effective_charge(z0, e1, t.z2);
                let f1: number = 0.3070722 * z1 * z1 * t.z2 / (b2 * a1 * t.a2);
                let f2: number = math.log(2.0 * this.ELECTRONMASS * b2 / t.iadj);
                if (this.switches.Shell) {
                    let etam2: number = 1.0 / (b * b * g * g);
                    let cadj: number = 1.0e-6 * t.iadj * t.iadj * etam2 * (0.422377 + etam2 * (0.0304043 - etam2 * 0.00038106)) + 1.0e-9 * t.iadj * t.iadj * t.iadj * etam2 * (3.858019 + etam2 * (-0.1667989 + etam2 * 0.00157955));
                    f2 -= cadj / t.z2;
                    if (this.switches.Leung) {
                        f2 -= (5.0 / 3.0) * math.log(2.0 * this.ELECTRONMASS * b2 / t.iadj) * (1.0e+03 * t.bind / (t.z2 * this.ELECTRONMASS)) - (t.iadj * t.iadj / (4.0 * this.ELECTRONMASS * this.ELECTRONMASS * b2));
                    }
                }
                let f6: number = 2.0 * math.log(g) - b2;
                // let f3: number = 0.0;
                let f3: number = this.lindhard(z1, a1, b);
                // console.log("f3 = " + f3.toString())
                let f4: number  = 1.0;
                if (this.switches.Barkas) {
                    let v: number = b * g / (this.ALPHA * math.sqrt(t.z2));
                    let fv: number = this.fva[9] * math.exp(-2.0 * math.log(v / 10.0));
                    f4 = 1.0 + 2.0 * z1 * fv / (math.sqrt(t.z2));
                }
                let f8: number = this.switches.Kinematic ? 0.5 * (-math.log(1.0 + 2.0 * (5.4858e-04 * g / a1)) - (5.4858e-04 * g / a1) * b2 / (g * g)) : 0.0;
                let f9: number = this.switches.Radiative ? (this.ALPHA / math.pi) * b2 * (6.0822 + math.log(2.0 * g) * (math.log(2.0 * g) * (2.4167 + 0.3333 * math.log(2.0 * g)) - 8.0314)) : 0.0;
                let Spa: number = 0.0;
                if ( this.switches.Pair ) {
                    let dpa: number = 1.0/math.sqrt(g);
                    let ldpa: number = math.log(dpa);
                    let l0: number = math.log(2.0*g);
                    // let Lpa0: number = (19.0/9.0)*(math.log(g/4.0) - 11.0/6.0);
                    let Lpa0s: number = (19.0/9.0)*math.log(183.0*math.exp(-1.0/3.0*math.log(t.z2))/(1.0 + 4.0*6.25470095193633*183.0*math.exp(-1.0/3.0*math.log(t.z2))/g));
                    let Lpa1: number = dpa*(4178.0/(81*math.pi*math.pi) - 21.0/27.0 - 248.0*l0/(27.0*math.pi*math.pi)
                        +(28.0*l0/9.0 - 446.0/27.0)*2.0*ldpa/(math.pi*math.pi) + 14.0*4.0*ldpa*ldpa/(9.0*math.pi*math.pi));
                    let Lpa: number = Lpa0s + Lpa1;
                    Spa=4.08803936906434e-06*(z1*z1/a1)*(t.z2*t.z2/t.a2)*(1.0 + 1.0/t.z2)*g*Lpa;
                }
                let Sbr: number = 0.0;
                if (this.switches.Bremsstrahlung ) {
                    let Bbr: number = math.log(1.0 + 2.0*g*0.179524783764566/(math.exp((1.0/3.0)*math.log(a1)) + math.exp((1.0/3.0)*math.log(t.a2)))/a1);
                    Sbr = 5.21721169334564e-07*(z1*z1/a1)*(z1*z1/a1)*(t.z2*t.z2/t.a2)*g*Bbr;
                }
                return f1 * (f2 * f4 + f3 + f6 - (delt / 2.0) + f8 + f9) + Sbr + Spa;
            }
        };
        let validateNumber = function(eventObject: JQuery.Event): boolean {
            if ($('#previous_target').val() !== 'Unknown') eventObject.data.first = false;
            let div: string = eventObject.data.name;
            if (eventObject.data.type === 'radio') {
                let name: string = div.split('_')[1];
                eventObject.data.valid = $("input[name=" + name + "]:checked").length === 1;
            } else {
                let input: JQuery<HTMLElement> = $('#' + div.split('_')[1].toUpperCase());
                let patt: RegExp = eventObject.data.type === 'int' ? new RegExp(/^[0-9]+$/i) : new RegExp(/^[0-9]+(\.[0-9]*|)(e[+-]?[0-9]+|)$/i);
                eventObject.data.valid = patt.test(input.val().toString());
                if (!eventObject.data.first) {
                    if (eventObject.data.valid) {
                        input.removeClass('is-invalid').addClass('is-valid');
                        $('#' + div + '_helpblock').html(eventObject.data.help);
                    } else {
                        input.removeClass('is-valid').addClass('is-invalid');
                        $('#' + div + '_helpblock').html(eventObject.data.help + " " + (eventObject.data.type === 'int' ? 'Integer' : 'Float') + " value required!");
                }
              }
            }
            eventObject.data.first = false;
            let validity: boolean[] = [];
            for (let j: number = 0; j < DeDx.validateDivs.length; j++) {
                validity.push(DeDx.validateDivs[j].valid);
            }
            let formvalid: boolean = validity.every(function(currentValue) { return currentValue; });
            $('#calculate').prop('disabled', !formvalid);
            return eventObject.data.valid;
        };
        let absorberTable = function(): void {
            let select: JQuery<HTMLElement> = $('#select_target');
            let previous_target: JQuery<HTMLElement> = $('#previous_target');
            let selected_target: string = previous_target.length === 1 ? previous_target.val().toString() : 'Unknown';
            let k: number = 0;
            for (let j: number = 0; j < DeDx.targets.length; j++) {
                if (DeDx.targets[j].name !== 'Unknown') {
                    let rowid: string = 't' + k;
                    k += 1;
                    $("<option id=\"" + rowid + "\"/>").html(DeDx.targets[j].name).appendTo(select);
                    if (DeDx.targets[j].name === selected_target) {
                        $("#" + rowid).prop('selected', true);
                    }
                }
            }
        };
        let calculate = function(eventObject: JQuery.Event): void {
            let task: string = $('input[name=task]:checked').val().toString();
            let re: number = Number($('#RE').val());
            let z0: number = Number($('#Z').val());
            let a1: number = Number($('#A').val());
            let target: string = $('#select_target').val().toString();
            // console.log("target = '" + target + "'");
            let bitmask: number = calculate_bitmask(eventObject);
            let type: string;
            let result: number;
            let unit: string;
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
            let r: JQuery<HTMLElement> = $('#result');
            let rr: JQuery<HTMLElement> = $("<tr id=\"rtr" + DeDx.nCalc + "\"/>");
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
        let calculate_bitmask = function(_eventObject: JQuery.Event): number {
            let bitmask: number = 0;
            let n: string;
            for (n in DeDx.switches) {
                if (!DeDx.switches.hasOwnProperty(n)) continue;
                let vv: boolean = $("#" + n).is(':checked');
                bitmask += vv ? Number($("#" + n).val()) : 0;
                DeDx.switches[n] = vv;
            }
            $('#bitmask').html(bitmask.toString());
            return bitmask
        };
        let resetForm = function(_eventObject: JQuery.Event): void {
            let recalc = document.getElementById('recalc') as HTMLFormElement;
            recalc.reset();
            $('#result').empty();
            DeDx.nCalc = 0;
            for (let j: number = 0; j < DeDx.validateDivs.length; j++) {
                if (DeDx.validateDivs[j].type !== "radio") {
                    let name: string = DeDx.validateDivs[j].name.split('_')[1].toUpperCase();
                    $("#" + name).removeClass('is-invalid').removeClass('is-valid');
                    $("#" + DeDx.validateDivs[j].name + "_helpblock").html(DeDx.validateDivs[j].help);
                }
            }
            let n: string;
            for (n in DeDx.switches) {
                if (!DeDx.switches.hasOwnProperty(n)) continue;
                DeDx.switches[n] = false
            }
            DeDx.switches.NewDelta = true;
            DeDx.switches.FiniteNuclearSize = true;
            $('#bitmask').html("40");
        };
        let csv = function(_eventObject: JQuery.Event): void {
            let rows: string[] = [];
            const header: string[] = ['ID', 'Task', 'E/R', 'Z', 'A', 'Target', 'Effects', 'Result', 'Units'];
            rows.push(header.join(','));
            let foo: JQuery<HTMLElement> = $('#result').children();
            if (foo.length === 0) {
                alert('No rows!');
                return;
            }
            for (let j: number = 0; j < foo.length; j++) {
                let r: string[] = [];
                for (let m: number = 0; m < foo[j].children.length; m++) {
                    r.push(foo[j].children[m].innerHTML.replace(/&nbsp;/g, ' '));
                }
                rows.push(r.join(','));
            }
            let c: string = rows.join('\r\n') + '\r\n';
            let download: JQuery<HTMLElement> = $("<a/>", {
                href: 'data:text/csv;charset=utf-8,' + encodeURIComponent(c),
                download: "crange.csv"
            });
            download.appendTo($("#downloadCSV"));
            download[0].click();
            $("#downloadCSV").empty();
        };
        if (DeDx.targets.length === 0) {
            $.getJSON("target.json", {}, function(data: Target[]) {
                return DeDx.targets = data;
            }).fail(function() {
                return alert("JSON error!");
            }).done(absorberTable);
        }
        for (let j: number = 0; j < DeDx.validateDivs.length; j++) {
            let div: Validation = DeDx.validateDivs[j];
            if ($("#" + div.name).length > 0) {
                $("#" + div.name).change(div, validateNumber).change();
            }
        }
        $('input[type=checkbox]').click(calculate_bitmask);
        $('#calculate').click(calculate);
        $('#resetForm').click(resetForm);
        $('#CSV').click(csv);
    }
);
