function togglePitch() {
        const type = document.querySelector('input[name="pitchType"]:checked').value;
        document.getElementById('identicalPitchInput').style.display = type === 'identical' ? 'block' : 'none';
        document.getElementById('unequalPitchInput').style.display = type === 'unequal' ? 'block' : 'none';
    }

    // Use provided Math.agm and Math.EllipticK/E implementations
/* The arithmetic-geometric mean of two non-negative numbers */
/*https://github.com/duetosymmetry/elliptic-integrals-js/blob/master/elliptic.js*/
Math.agm = function(a0,g0) {
  var maxIter = 50;
  var an = (a0+g0)/2;
  var gn = Math.sqrt(a0*g0);
  var iter;
  for (iter=0; iter<maxIter && Math.abs(an-gn)>1e-15; iter++) {
    a0 = 0.5*(an+gn);
    g0 = Math.sqrt(an*gn);
    an = a0;
    gn = g0;
  }
  if (iter==maxIter) console.warn("Math.agm hit iteration limit, may not have converged");
  return an;
};

/* Complete elliptic integral of the first kind K(m) */
Math.EllipticK = function(m) {
  var kprime = Math.sqrt(1-m);
  return 0.5*Math.PI/Math.agm(1, kprime);
};

/* Complete elliptic integral of the second kind E(m) */
Math.EllipticE = function(m) {
  var maxIter=50, iter=0;
  var kprime = Math.sqrt(1-m);
  var a0=1, g0=kprime, an=a0, gn=g0;
  var twoPow=0.25;
  var partialSum = 1 - 0.5*m;
  while(Math.abs(an-gn)>1e-15 && iter<maxIter) {
    partialSum -= twoPow*(an-gn)*(an-gn);
    twoPow *= 2;
    a0 = 0.5*(an+gn);
    g0 = Math.sqrt(an*gn);
    an = a0; gn = g0;
    iter++;
  }
  if (iter==maxIter) console.warn("Math.EllipticE hit iteration limit, may not have converged");
  return 0.5*Math.PI*partialSum/an;
};

function f_L_self_helical(r0, N, rout, p) {
        const mu0 = 4 * Math.PI * 1e-7;
        // build z array length N
        const z = [0];
        for (let i = 1; i < N; i++) {
            z[i] = z[i-1] + p[i-1];
        }
        const M = Array.from({length: N}, () => Array(N).fill(0));
        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                if (i === j) {
                    M[i][j] = mu0 * rout * (Math.log(8 * rout / r0) - 1.75);
                } else {
                    const dz = z[i] - z[j];
                    let m = 4 * rout * rout / (4 * rout * rout + dz * dz);
                    m = Math.min(Math.max(m, 1e-12), 1 - 1e-12);
                    const K = Math.EllipticK(m);
                    const E = Math.EllipticE(m);
                    M[i][j] = mu0 * Math.sqrt(rout * rout / m) * ((2 - m) * K - 2 * E);
                }
            }
        }
        let Lsum = 0;
        for (let ii = 0; ii < N; ii++) {
            for (let jj = 0; jj < N; jj++) {
                Lsum += M[ii][jj];
            }
        }
        return Lsum;
    }

    function f_L_self_spiral(r0, N, rout, p) {
        const mu0 = 4 * Math.PI * 1e-7;
        // build r array length N
        const r = [rout];
        for (let i = 1; i < N; i++) {
            r[i] = rout - p.slice(0, i).reduce((s, x) => s + x, 0);
        }
        const M = Array.from({length: N}, () => Array(N).fill(0));
        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                if (i === j) {
                    M[i][j] = mu0 * r[i] * (Math.log(8 * r[i] / r0) - 1.75);
                } else {
                    let m = 4 * r[i] * r[j] / ((r[i] + r[j]) * (r[i] + r[j]));
                    m = Math.min(Math.max(m, 1e-12), 1 - 1e-12);
                    const K = Math.EllipticK(m);
                    const E = Math.EllipticE(m);
                    M[i][j] = mu0 * Math.sqrt(r[i] * r[j] / m) * ((2 - m) * K - 2 * E);
                }
            }
        }
        let Lsum = 0;
        for (let ii = 0; ii < N; ii++) {
            for (let jj = 0; jj < N; jj++) {
                Lsum += M[ii][jj];
            }
        }
        return Lsum;
    }

 // Symmetric magnetic field component
    function f_H_sym(R, p, m, N) {
        const xx = N % 2;
        const NN = xx === 0 ? N/2 : (N+1)/2;
        if (m === 1 || m === N) return 0;
        let pp_left = [], pp_right = [];
        if (m <= NN) {
            for (let k = m-1, sumL=0; k >= 1; k--) { sumL += p[k-1]; pp_left.push(sumL); }
            for (let k = m, sumR=0; k < m + pp_left.length; k++) { sumR += p[k-1]; pp_right.push(sumR); }
        } else {
            const lower = 2*m - N;
            for (let k = m-1, sumL=0; k >= lower; k--) { sumL += p[k-1]; pp_left.push(sumL); }
            for (let k = m, sumR=0; k < m + pp_left.length; k++) { sumR += p[k-1]; pp_right.push(sumR); }
        }
        return pp_left.reduce((Hsum, _, i) => {
            const x1 = pp_left[i]/R, x2 = pp_right[i]/R;
            const term = Math.sqrt((1 + x1*x1)/Math.pow(x1*x1 - 1, 2)
                + (1 + x2*x2)/Math.pow(x2*x2 - 1, 2)
                - 2*(x1*x2 - 1)/((x1*x1 - 1)*(x2*x2 - 1)))
                / (2 * Math.PI * R);
            return Hsum + term;
        }, 0);
    }

    // Asymmetric magnetic field component
    function f_H_asy(R, p, m, N) {
        const xx = N % 2;
        const NN = xx === 0 ? N/2 : (N+1)/2;
        if (m === NN && xx === 1) return 0;
        let Hsum = 0;
        if (m <= NN) {
            let offset = 0;
            for (let k = m; k < m + (m-1); k++) offset += p[k-1];
            let pp_right = [];
            for (let k = 2*m - 1, sumR=0; k <= N-1; k++) { sumR += p[k-1]; pp_right.push(sumR + offset); }
            Hsum = pp_right.reduce((sum,val)=> sum + val/(val*val + R*R)/(2*Math.PI), 0);
        } else if (m < N) {
            let offset = 0;
            for (let k = m-1; k >= 2*m - N; k--) offset += p[k-1];
            let pp_left = [];
            for (let k = 2*m - N - 1, sumL=0; k >= 1; k--) { sumL += p[k-1]; pp_left.push(sumL + offset); }
            Hsum = pp_left.reduce((sum,val)=> sum + val/(val*val + R*R)/(2*Math.PI), 0);
        } else {
            let sumL=0, pp_left=[];
            for (let k = 2*m - N - 1; k >= 1; k--) { sumL += p[k-1]; pp_left.push(sumL); }
            Hsum = pp_left.reduce((sum,val)=> sum + val/(val*val + R*R)/(2*Math.PI), 0);
        }
        return Hsum;
    }



    function calculate() {
        const OD_mm = parseFloat(document.getElementById('OD').value);
        const wire_r_mm = parseFloat(document.getElementById('r0').value);
        const coil_r_max_mm = OD_mm/2;
        const shape = document.querySelector('input[name="coilShape"]:checked').value;
        const f = parseFloat(document.getElementById('f').value) * 1e6;
        const N = parseInt(document.getElementById('N').value);
        const sigma = parseFloat(document.getElementById('s').value) * 1e7;

        let p_mm = [];
        const type = document.querySelector('input[name="pitchType"]:checked').value;
        if (type === 'identical') {
            const p0 = parseFloat(document.getElementById('p0').value);
            p_mm = Array(N-1).fill(p0);
        } else {
            p_mm = document.getElementById('pList').value.split(',').map(v=>parseFloat(v));
        }

        let coilLength_mm;
        if (shape === 'spiral') {
            const radii = [coil_r_max_mm];
            for (let i = 1; i < N; i++) radii.push(radii[i-1] - p_mm[i-1]);
            coilLength_mm = radii.reduce((sum, r) => sum + 2*Math.PI*r, 0);
        } else {
            coilLength_mm = Math.PI * OD_mm * N; // helical length = π·OD·N
        }

        const coilLength_m = coilLength_mm/1000;
        const area = Math.PI * Math.pow(wire_r_mm/1000, 2);
        const DCR = coilLength_m / (area * sigma);

        const r0 = wire_r_mm/1000;
        const p_list_m = p_mm.map(v=>v/1000);
        const Hm = [];
        for (let m=1; m<=N; m++) {
            Hm.push(f_H_asy(r0, p_list_m, m, N) + f_H_sym(r0, p_list_m, m, N));
        }
        const mu0 = 4 * Math.PI * 1e-7;
        const sd = 1 / Math.sqrt(Math.PI * f * mu0 * sigma);
        const x = 2 * r0 / sd;
        const sumH2 = Hm.reduce((s,h)=>s+h*h,0);
        const Gp = (r0 - sd > 0)
            ? (8 * sumH2 * Math.PI**2 * sd**2 * x**3 * (x-1) / (((2*x+1)**2 + 2))) / N
            : (Math.PI**2 * r0**2 * sumH2 * (x/2)**4 / (1 + (x/2)**4/48)) / N;
        const Rskin = (r0 - sd > 0)
            ? (0.25 + r0/(2*sd) + 3*sd/(32*r0)) / (Math.PI * r0**2 * sigma)*coilLength_m
            : (1 + Math.pow(r0/sd,4)/48) / (Math.PI * r0**2 * sigma)*coilLength_m;
        const Rprox = Gp * Rskin;
        const Rtotal = (1 + Gp) * Rskin;

        document.getElementById('coilLength').textContent = coilLength_mm.toFixed(3);
        document.getElementById('DCR').textContent = DCR.toFixed(6);
        document.getElementById('Rskin').textContent = Rskin.toFixed(6);
        document.getElementById('Rprox').textContent = Rprox.toFixed(6);
        document.getElementById('Rtotal').textContent = Rtotal.toFixed(6);
        document.getElementById('Gp').textContent = Gp.toFixed(6);


        // Inductance
        const r0_m = wire_r_mm/1000;
        const rout_m = coil_r_max_mm/1000;
        const p_m = p_mm.map(v=>v/1000);
        let L_self = (shape === 'helical')
            ? f_L_self_helical(r0_m, N, rout_m, p_m)
            : f_L_self_spiral(r0_m, N, rout_m, p_m);
        const L_uH = L_self * 1e6;

        document.getElementById('L').textContent = L_uH.toFixed(3);
        // Other outputs updated similarly...
    }
