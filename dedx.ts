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

interface RangeTable {
    z1: number;
    a1: number;
    bitmask: number;
    target: string;
    range: number[];
    interpolateRange(energy: number): number;
    interpolateEnergy(energy: number, r0: number): number;
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
    rangeTables: RangeTable[];
    absorber(target: string): Target;
    effective_charge(z0: number, e1: number, z2: number): number;
    delta(g: number, t: Target): number;
    olddelta(g: number, t: Target): number;
    lngamma(z: math.Complex): math.Complex;
    hyperg(a: math.Complex, b: math.Complex, z: math.Complex): math.Complex;
    lindhard(zz: number, aa: number, bb: number): number;
    stop(e1: number, z0: number, a1: number, target: string): number;
    range(e: number, z1: number, a1: number, bitmask: number, target: string): number;
    benton(e: number, z1: number, a1: number, target: string): number;
    integrate(i: number, z1: number, a1: number, target: string): number;
    renergy(e: number, r0: number, z1: number, a1: number, bitmask: number, target: string): number;
}

declare namespace JQuery {
    interface Event {
        data: Validation;
    }
}

$(
    function(): void {
        const MAXE: number = 200;
        const LOGTENEMIN: number = 0;
        const LOGTENEMAX: number = 6;
        const M_LN10: number = 2.30258509299404568402;
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
                name: "div_task",
                type: "radio",
                valid: false,
                first: true,
                help: ""
            }, {
                name: "div_re",
                type: "float",
                valid: false,
                first: true,
                help: "Units are A&nbsp;MeV."
            }, {
                name: "div_dx",
                type: "float",
                valid: false,
                first: true,
                help: "Distance traversed in target material. Units are g&nbsp;cm<sup>-2</sup>."
            }, {
                name: "div_z",
                type: "int",
                valid: false,
                first: true,
                help: "Number of protons; atomic number."
            }, {
                name: "div_a",
                type: "int",
                valid: false,
                first: true,
                help: "Number of nucleons; atomic mass."
            }],
            rangeTables: [],
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
            },
            range: function(e: number, z1: number, a1: number, bitmask: number, target: string): number {
                // const t: Target = this.absorber(target);
                let rt: RangeTable;
                for (let i: number = 0; i < this.rangeTables.length; i++) {
                    rt = this.rangeTables[i];
                    if ((z1 == rt.z1) && (a1 == rt.a1) && (bitmask == rt.bitmask) && (target == rt.target)) {
                        return rt.interpolateRange(e);
                    }
                }
                rt = CreateRangeTable(z1, a1, bitmask, target);
                this.rangeTables.push(rt);
                return rt.interpolateRange(e);
            },
            benton: function(e: number, z1: number, a1: number, target: string): number {
                const t: Target = this.absorber(target);
                const amn: number[][][] = [[[-8.72500,  1.88000,  0.741900,  0.752000],
                                            [ 0.83090,  0.11140, -0.528800, -0.555890],
                                            [-0.13396, -0.06481,  0.126423,  0.128431],
                                            [ 0.01262,  0.00540, -0.009341, -0.009306]],
                                           [[-7.6604e-01,  2.5398e+00, -2.4598e-01,  0.0000e+00],
                                            [ 7.3736e-02, -3.1200e-01,  1.1548e-01,  0.0000e+00],
                                            [ 4.0556e-02,  1.8664e-02, -9.9661e-03,  0.0000e+00],
                                            [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00]],
                                           [[-8.0155e+00,  1.8371e+00,  4.5233e-02, -5.9898e-03],
                                            [ 3.6916e-01, -1.4520e-02, -9.5873e-04, -5.2315e-04],
                                            [-1.4307e-02, -3.0142e-02,  7.1303e-03, -3.3802e-04],
                                            [ 3.4718e-03,  2.3603e-03, -6.8538e-04,  3.9405e-05]]];
                const czmn: number[][] = [[-6.000e-05,  5.252e-02,  1.285e-01,  0.000e+00],
                                          [-1.850e-03,  7.355e-02,  7.171e-02, -2.723e-02],
                                          [-7.930e-02,  3.323e-01, -1.234e-01,  1.530e-02],
                                          [ 2.200e-01,  0.000e+00,  0.000e+00,  0.000e+00]];
                const cjoin: number[][] = [[ 0.94,  20.19, -84.08,  132.98, -30.77, -102.29, 64.03],
                                           [12.62, -51.96, 199.14, -367.09, 327.06, -108.57,  0.00]];
                const cr: number = this.PROTONMASS/this.ATOMICMASSUNIT;
                //
                // Compute join[4].
                //
                let join: number[] = [0, 0, 0, 0];
                for (let l: number = 0; l < 2; l++) {
                    let m: number = 6;
                    join[l] = cjoin[l][m];
                    while (m > 0) join[l] = 0.001*t.iadj*join[l] + cjoin[l][--m];
                }
                //
                // Compute prnglo[3].
                //
                const logt: number = math.log( e * cr );
                const logi: number = math.log( t.iadj );
                let prnglo: number[] = [0, 0, 0];
                for (let l: number = 0; l < 3; l++) {
                    let loglambda: number = 0.0;
                    for (let m: number = 3; m >= 0; m--) {
                        let n: number = 3;
                        let term: number = amn[l][m][n];
                        while (n > 0) term = term*logt + amn[l][m][--n];
                        loglambda += term;
                        loglambda *= logi;
                    }
                    loglambda /= logi;
                    loglambda += math.log(t.a2/t.z2);
                    prnglo[l] = math.exp(loglambda);
                    if (l == 1) prnglo[l] *= 1.0e-03;
                }
                //
                // Compute cz[4].
                //
                const g: number = 1.0 + e/this.ATOMICMASSUNIT;
                const b: number = math.sqrt(1.0 -1.0/(g*g));
                const x: number = 137.0*b/z1;
                let cz: number[] = [0, 0, 0, 0];
                for (let m: number = 0; m < 4; m++) {
                    cz[m]=0.0;
                    for (let n: number = 3; n >= 0; n--) {
                        cz[m] += czmn[m][n];
                        cz[m] *= x;
                    }
                    cz[m] /= x;
                }
                let l: number = 1;
                if (e < join[0]) l = 0;
                if (e > join[1]) l = 2;
                let n: number;
                if (x <= 0.2) {
                    n = 0;
                } else if ( x > 0.2 && x <= 2.0 ) {
                    n = 1;
                } else if ( x > 2.0 && x <= 3.0 ) {
                    n = 2;
                } else { // if ( x > 3.0 ) {
                    n = 3;
                }
                let bzz: number = (31.8 + 3.86*math.exp((5.0/8.0)*logi)) * (t.a2/t.z2)*1.0e-06*math.exp((8.0/3.0)*math.log(z1));
                return (( (a1/cr)/(z1*z1) )*(prnglo[l] + bzz*cz[n]));
            },
            integrate: function(i: number, z1: number, a1: number, target: string): number {
                // const t: Target = this.absorber(target);
                // const rel: number = 0.0;
                const e0: number = energyTable(i-1);
                const de2: number = (energyTable(i) - e0)/2.0;
                const e1: number = e0 + 1.33998104*de2;
                const dedx1: number = this.stop(e1, z1, a1, target);
                const e2: number = e0 + 1.86113631*de2;
                const dedx2: number = this.stop(e2, z1, a1, target);
                const e3: number = e0 + 0.13886369*de2;
                const dedx3: number = this.stop(e3, z1, a1, target);
                const e4: number = e0 + 0.66001869*de2;
                const dedx4: number = this.stop(e4, z1, a1, target);
                return de2 * (0.65214515/dedx1 +
                              0.34785485/dedx2 +
                              0.34785485/dedx3 +
                              0.65214515/dedx4);
            },
            renergy: function(e: number, r0: number, z1: number, a1: number, bitmask: number, target: string): number {
                let rt: RangeTable;
                for (let i: number = 0; i < this.rangeTables.length; i++) {
                    rt = this.rangeTables[i];
                    if ((z1 == rt.z1) && (a1 == rt.a1) && (bitmask == rt.bitmask) && (target == rt.target)) {
                        return rt.interpolateEnergy(e, r0);
                    }
                }
                rt = CreateRangeTable(z1, a1, bitmask, target);
                this.rangeTables.push(rt);
                return rt.interpolateEnergy(e, r0);
            }
        };
        let CreateRangeTable = function(z1: number, a1: number, bitmask: number, target: string): RangeTable {
            let i: number = 0;
            let range: number[] = [];
            while (energyTable(i) < 8.0) {
                range.push(DeDx.benton(energyTable(i), z1, a1, target));
                i++;
            }
            while (i < MAXE) {
                range.push(range[i-1] + DeDx.integrate(i, z1, a1, target));
                i++;
            }
            let rt: RangeTable = {
                z1: z1,
                a1: a1,
                bitmask: bitmask,
                target: target,
                range: range,
                interpolateRange: function(energy: number): number {
                    if (energy > energyTable(0)) {
                        let i: number = 1;
                        while (energy > energyTable(i)) i++;
                        return (this.range[i-1] + (energy - energyTable(i-1)) *
                                (this.range[i]-this.range[i-1])/(energyTable(i)-energyTable(i-1)));
                    } else {
                        return energy*this.range[0]/energyTable(0);
                    }
                },
                interpolateEnergy: function(energy: number, r0: number): number {
                    let r: number;
                    // Not sure why energy would ever be negative.
                    if ( energy > 0.0 ){
                        //
                        // Compute the expected total range at this energy.
                        //
                        let rr: number = this.interpolateRange(energy);
                        r = rr - r0;
                        if (r < 0.0) return(0.0);
                    } else {
                        r = r0;
                    }
                    if (r > range[0]) {
                        let i: number = 1;
                        while (r > this.range[i]) i++;
                        return (energyTable(i-1) + (r - this.range[i-1]) * (energyTable(i)-energyTable(i-1)) / (this.range[i]-this.range[i-1]));
                    } else {
                        return energyTable(0) * r / this.range[0];
                    }
                }
            };
            return rt;
        };
        let energyTable = function(i: number): number {
            return math.exp(M_LN10*(LOGTENEMIN + (i)*(LOGTENEMAX-LOGTENEMIN)/(MAXE - 1.0)));
        };
        let validateNumber = function(eventObject: JQuery.Event): boolean {
            if ($("#previous_target").val() !== "Unknown") eventObject.data.first = false;
            let div: string = eventObject.data.name;
            if (eventObject.data.type === "radio") {
                let name: string = div.split("_")[1];
                eventObject.data.valid = $("input[name=" + name + "]:checked").length === 1;
            } else {
                let input: JQuery<HTMLElement> = $("#" + div.split("_")[1].toUpperCase());
                let patt: RegExp = eventObject.data.type === "int" ? new RegExp(/^[0-9]+$/i) : new RegExp(/^[0-9]+(\.[0-9]*|)(e[+-]?[0-9]+|)$/i);
                eventObject.data.valid = patt.test(input.val().toString());
                if (!eventObject.data.first) {
                    if (eventObject.data.valid) {
                        input.removeClass("is-invalid").addClass("is-valid");
                        $("#" + div + "_helpblock").html(eventObject.data.help);
                    } else {
                        input.removeClass("is-valid").addClass("is-invalid");
                        $("#" + div + "_helpblock").html(eventObject.data.help + " " + (eventObject.data.type === "int" ? "Integer" : "Float") + " value required!");
                }
              }
            }
            eventObject.data.first = false;
            let validity: boolean[] = [];
            for (let j: number = 0; j < DeDx.validateDivs.length; j++) {
                validity.push(DeDx.validateDivs[j].valid);
            }
            let formvalid: boolean = validity.every(function(currentValue) { return currentValue; });
            // $("#calculate").prop("disabled", !formvalid);
            if (formvalid) {
                $("#calculate").removeAttr("disabled");
            } else {
                $("#calculate").attr("disabled", "disabled");
            }
            return eventObject.data.valid;
        };
        let absorberTable = function(): void {
            let select: JQuery<HTMLElement> = $("#select_target");
            let previous_target: JQuery<HTMLElement> = $("#previous_target");
            let selected_target: string = previous_target.length === 1 ? previous_target.val().toString() : "Unknown";
            let k: number = 0;
            for (let j: number = 0; j < DeDx.targets.length; j++) {
                if (DeDx.targets[j].name !== "Unknown") {
                    let rowid: string = "t" + k;
                    k += 1;
                    $("<option id=\"" + rowid + "\"/>").html(DeDx.targets[j].name).appendTo(select);
                    if (DeDx.targets[j].name === selected_target) {
                        $("#" + rowid).attr("selected", "selected");
                    }
                }
            }
        };
        let task_help = function(_eventObject: JQuery.Event): void {
            let task: string = $("input[name=task]:checked").val().toString();
            if (task == "e"){
                $("#div_dx").removeClass("d-none").addClass("d-block");
                // $("#DX").attr("value", "");
                $("#DX").val("");
            } else {
                $("#div_dx").removeClass("d-block").addClass("d-none");
                // $("#DX").attr("value", "0");
                $("#DX").val("0");
            }
        };
        let calculate = function(eventObject: JQuery.Event): void {
            $("#calculate").attr("disabled", "disabled");
            let task: string = $("input[name=task]:checked").val().toString();
            let re: number = Number($("#RE").val());
            let dx: number = Number($("#DX").val());
            let z0: number = Number($("#Z").val());
            let a1: number = Number($("#A").val());
            let target: string = $("#select_target").val().toString();
            let bitmask: number = calculate_bitmask(eventObject);
            let type: string;
            let result: number;
            let unit: string;
            switch (task) {
                case "r":
                    type = "Range";
                    console.log("DeDx.range(" + re + ", " + z0 + ", " + a1 + ", " + bitmask + ", '" + target + "');");
                    result = DeDx.range(re, z0, a1, bitmask, target);
                    unit = "g&nbsp;cm<sup>-2</sup>";
                    break;
                case "e":
                    type = "Energy";
                    console.log("DeDx.renergy(" + re + ", " + dx + ", " + z0 + ", " + a1 + ", " + bitmask + ", '" + target + "');");
                    result = DeDx.renergy(re, dx, z0, a1, bitmask, target);
                    unit = "A&nbsp;MeV";
                    break;
                case "d":
                    type = "d<var>E</var>/d<var>x</var>";
                    console.log("DeDx.stop(" + re + ", " + z0 + ", " + a1 + ", '" + target + "');");
                    result = DeDx.stop(re, z0, a1, target);
                    unit = "A&nbsp;MeV&nbsp;g<sup>-1</sup>&nbsp;cm<sup>2</sup>";
            }
            DeDx.nCalc += 1;
            let r: JQuery<HTMLElement> = $("#result");
            let rr: JQuery<HTMLElement> = $("<tr id=\"rtr" + DeDx.nCalc + "\"/>");
            $("<td/>").html(DeDx.nCalc.toString()).appendTo(rr);
            $("<td/>").html(type).appendTo(rr);
            $("<td/>").html(re.toString()).appendTo(rr);
            $("<td/>").html(dx.toString()).appendTo(rr);
            $("<td/>").html(z0.toString()).appendTo(rr);
            $("<td/>").html(a1.toString()).appendTo(rr);
            $("<td/>").html(target).appendTo(rr);
            $("<td/>").html(bitmask.toString()).appendTo(rr);
            $("<td/>").html(result.toString()).appendTo(rr);
            $("<td/>").html(unit).appendTo(rr);
            rr.appendTo(r);
            $("#calculate").removeAttr("disabled").html("Calculate");
        };
        let calculate_bitmask = function(_eventObject: JQuery.Event): number {
            let bitmask: number = 0;
            let n: string;
            for (n in DeDx.switches) {
                if (!DeDx.switches.hasOwnProperty(n)) continue;
                let vv: boolean = $("#" + n).is(":checked");
                bitmask += vv ? Number($("#" + n).val()) : 0;
                DeDx.switches[n] = vv;
            }
            $("#bitmask").html(bitmask.toString());
            return bitmask
        };
        let resetForm = function(_eventObject: JQuery.Event): void {
            let recalc = document.getElementById("recalc") as HTMLFormElement;
            recalc.reset();
            $("#result").empty();
            DeDx.nCalc = 0;
            for (let j: number = 0; j < DeDx.validateDivs.length; j++) {
                if (DeDx.validateDivs[j].type !== "radio") {
                    let name: string = DeDx.validateDivs[j].name.split("_")[1].toUpperCase();
                    $("#" + name).removeClass("is-invalid").removeClass("is-valid");
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
            $("#bitmask").html("40");
        };
        let csv = function(_eventObject: JQuery.Event): void {
            let rows: string[] = [];
            const header: string[] = ["ID", "Task", "E", "dx", "Z", "A", "Target", "Effects", "Result", "Units"];
            rows.push(header.join(","));
            let foo: JQuery<HTMLElement> = $("#result").children();
            if (foo.length === 0) {
                alert("No rows!");
                return;
            }
            for (let j: number = 0; j < foo.length; j++) {
                let r: string[] = [];
                for (let m: number = 0; m < foo[j].children.length; m++) {
                    r.push(foo[j].children[m].innerHTML.replace(/&nbsp;/g, " "));
                }
                rows.push(r.join(","));
            }
            let c: string = rows.join("\r\n") + "\r\n";
            let download: JQuery<HTMLElement> = $("<a/>", {
                href: "data:text/csv;charset=utf-8," + encodeURIComponent(c),
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
        $("input[type=radio]").on("click", task_help);
        $("input[type=checkbox]").on("click", calculate_bitmask);
        $("#calculate").on("mousedown", function(): void {
            let spin: string = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>Calculating...';
            $("#calculate").html(spin);
        }).on("mouseup", calculate);
        $("#resetForm").on("click", resetForm);
        $("#CSV").on("click", csv);
    }
);
