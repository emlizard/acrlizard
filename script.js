document.addEventListener('DOMContentLoaded', () => {
    const App = {
        elements: {
            themeToggle: document.getElementById('theme-toggle'),
            calculateBtn: document.getElementById('calculateBtn'),
            pitchTypeRadios: document.querySelectorAll('input[name="pitchType"]'),
            identicalPitchInput: document.getElementById('identicalPitchInput'),
            unequalPitchInput: document.getElementById('unequalPitchInput'),
            // Input fields
            f: document.getElementById('f'),
            s: document.getElementById('s'),
            r0: document.getElementById('r0'),
            N: document.getElementById('N'),
            OD: document.getElementById('OD'),
            p0: document.getElementById('p0'),
            pList: document.getElementById('pList'),
            coilShapeRadios: document.querySelectorAll('input[name="coilShape"]'),
            // Output fields
            coilLength: document.getElementById('coilLength'),
            DCR: document.getElementById('DCR'),
            Rskin: document.getElementById('Rskin'),
            Rprox: document.getElementById('Rprox'),
            Rtotal: document.getElementById('Rtotal'),
            L: document.getElementById('L'),
            Gp: document.getElementById('Gp'),
        },

        init() {
            this.initTheme();
            this.initEventListeners();
            this.togglePitch(); // Set initial state
            this.calculate();   // Initial calculation on load
        },

        initTheme() {
            const sunIcon = `☀️`;
            const moonIcon = `🌙`;
            const setTheme = (theme) => {
                document.body.setAttribute('data-theme', theme);
                this.elements.themeToggle.innerHTML = theme === 'dark' ? sunIcon : moonIcon;
                localStorage.setItem('theme', theme);
            };
            this.elements.themeToggle.addEventListener('click', () => {
                const currentTheme = document.body.getAttribute('data-theme');
                setTheme(currentTheme === 'light' ? 'dark' : 'light');
            });
            const savedTheme = localStorage.getItem('theme') || 'light';
            setTheme(savedTheme);
        },

        initEventListeners() {
            this.elements.calculateBtn.addEventListener('click', this.calculate.bind(this));
            this.elements.pitchTypeRadios.forEach(radio => {
                radio.addEventListener('change', this.togglePitch.bind(this));
            });
        },

        togglePitch() {
            const type = document.querySelector('input[name="pitchType"]:checked').value;
            this.elements.identicalPitchInput.style.display = type === 'identical' ? 'block' : 'none';
            this.elements.unequalPitchInput.style.display = type === 'unequal' ? 'block' : 'none';
        },

        // --- CORE CALCULATION LOGIC (PRESERVED FROM ORIGINAL) ---
        agm(a0, g0) {
            let maxIter = 50, an = (a0 + g0) / 2, gn = Math.sqrt(a0 * g0), iter;
            for (iter = 0; iter < maxIter && Math.abs(an - gn) > 1e-15; iter++) {
                a0 = 0.5 * (an + gn); g0 = Math.sqrt(an * gn); an = a0; gn = g0;
            }
            if (iter === maxIter) console.warn("Math.agm hit iteration limit, may not have converged");
            return an;
        },
        EllipticK(m) {
            if (m >= 1) return Infinity;
            const kprime = Math.sqrt(1 - m);
            return 0.5 * Math.PI / this.agm(1, kprime);
        },
        EllipticE(m) {
            if (m > 1) return 0;
            let maxIter = 50, iter = 0, a0 = 1, g0 = Math.sqrt(1 - m), an = a0, gn = g0, twoPow = 0.25;
            let partialSum = 1 - 0.5 * m;
            while (Math.abs(an - gn) > 1e-15 && iter < maxIter) {
                partialSum -= twoPow * (an - gn) * (an - gn);
                twoPow *= 2; a0 = 0.5 * (an + gn); g0 = Math.sqrt(an * gn); an = a0; gn = g0; iter++;
            }
            if (iter === maxIter) console.warn("Math.EllipticE hit iteration limit, may not have converged");
            return 0.5 * Math.PI * partialSum / an;
        },
        f_L_self_helical(r0, N, rout, p) {
            const mu0 = 4 * Math.PI * 1e-7;
            const z = [0];
            for (let i = 1; i < N; i++) z[i] = z[i - 1] + p[i - 1];
            let Lsum = 0;
            for (let i = 0; i < N; i++) {
                for (let j = 0; j < N; j++) {
                    if (i === j) {
                        Lsum += mu0 * rout * (Math.log(8 * rout / r0) - 1.75);
                    } else {
                        const dz = z[i] - z[j];
                        let m = 4 * rout * rout / (4 * rout * rout + dz * dz);
                        m = Math.min(Math.max(m, 1e-12), 1 - 1e-12);
                        const K = this.EllipticK(m);
                        const E = this.EllipticE(m);
                        Lsum += mu0 * Math.sqrt(rout * rout / m) * ((2 - m) * K - 2 * E);
                    }
                }
            }
            return Lsum;
        },
        f_L_self_spiral(r0, N, rout, p) {
            const mu0 = 4 * Math.PI * 1e-7;
            const r = [rout];
            for (let i = 1; i < N; i++) r[i] = rout - p.slice(0, i).reduce((s, x) => s + x, 0);
            let Lsum = 0;
            for (let i = 0; i < N; i++) {
                for (let j = 0; j < N; j++) {
                    if (i === j) {
                        Lsum += mu0 * r[i] * (Math.log(8 * r[i] / r0) - 1.75);
                    } else {
                        let m = 4 * r[i] * r[j] / ((r[i] + r[j]) ** 2);
                        m = Math.min(Math.max(m, 1e-12), 1 - 1e-12);
                        const K = this.EllipticK(m);
                        const E = this.EllipticE(m);
                        Lsum += mu0 * Math.sqrt(r[i] * r[j] / m) * ((2 - m) * K - 2 * E);
                    }
                }
            }
            return Lsum;
        },
        f_H_sym(R, p, m, N) {
            const xx = N % 2; const NN = xx === 0 ? N / 2 : (N + 1) / 2;
            if (m === 1 || m === N) return 0;
            let pp_left = [], pp_right = [];
            if (m <= NN) {
                for (let k = m - 1, sumL = 0; k >= 1; k--) { sumL += p[k - 1]; pp_left.push(sumL); }
                for (let k = m, sumR = 0; k < m + pp_left.length; k++) { sumR += p[k - 1]; pp_right.push(sumR); }
            } else {
                const lower = 2 * m - N;
                for (let k = m - 1, sumL = 0; k >= lower; k--) { sumL += p[k - 1]; pp_left.push(sumL); }
                for (let k = m, sumR = 0; k < m + pp_left.length; k++) { sumR += p[k - 1]; pp_right.push(sumR); }
            }
            return pp_left.reduce((Hsum, _, i) => {
                const x1 = pp_left[i] / R, x2 = pp_right[i] / R;
                const term = Math.sqrt((1 + x1 * x1) / ((x1 * x1 - 1) ** 2) + (1 + x2 * x2) / ((x2 * x2 - 1) ** 2) - 2 * (x1 * x2 - 1) / ((x1 * x1 - 1) * (x2 * x2 - 1))) / (2 * Math.PI * R);
                return Hsum + term;
            }, 0);
        },
        f_H_asy(R, p, m, N) {
            const xx = N % 2; const NN = xx === 0 ? N / 2 : (N + 1) / 2;
            if (m === NN && xx === 1) return 0;
            let Hsum = 0;
            if (m <= NN) {
                let offset = 0;
                for (let k = m; k < m + (m - 1); k++) offset += p[k - 1];
                let pp_right = [];
                for (let k = 2 * m - 1, sumR = 0; k <= N - 1; k++) { sumR += p[k - 1]; pp_right.push(sumR + offset); }
                Hsum = pp_right.reduce((sum, val) => sum + val / (val * val + R * R) / (2 * Math.PI), 0);
            } else if (m < N) {
                let offset = 0;
                for (let k = m - 1; k >= 2 * m - N; k--) offset += p[k - 1];
                let pp_left = [];
                for (let k = 2 * m - N - 1, sumL = 0; k >= 1; k--) { sumL += p[k - 1]; pp_left.push(sumL + offset); }
                Hsum = pp_left.reduce((sum, val) => sum + val / (val * val + R * R) / (2 * Math.PI), 0);
            } else {
                let sumL = 0, pp_left = [];
                for (let k = 2 * m - N - 1; k >= 1; k--) { sumL += p[k - 1]; pp_left.push(sumL); }
                Hsum = pp_left.reduce((sum, val) => sum + val / (val * val + R * R) / (2 * Math.PI), 0);
            }
            return Hsum;
        },

        calculate() {
            try {
                // Get all values
                const OD_mm = parseFloat(this.elements.OD.value);
                const wire_r_mm = parseFloat(this.elements.r0.value);
                const coil_r_max_mm = OD_mm / 2;
                const shape = document.querySelector('input[name="coilShape"]:checked').value;
                const f = parseFloat(this.elements.f.value) * 1e6;
                const N = parseInt(this.elements.N.value);
                const sigma = parseFloat(this.elements.s.value) * 1e7;

                let p_mm = [];
                const type = document.querySelector('input[name="pitchType"]:checked').value;
                if (type === 'identical') {
                    p_mm = Array(N - 1).fill(parseFloat(this.elements.p0.value));
                } else {
                    p_mm = this.elements.pList.value.split(',').map(v => parseFloat(v.trim()));
                }

                if (p_mm.some(isNaN)) throw new Error("Pitches must be valid numbers.");
                if (p_mm.length !== N - 1 && N > 1) throw new Error(`Number of pitches (${p_mm.length}) must be N-1 (${N - 1}).`);

                // Calculations
                let coilLength_mm;
                if (shape === 'spiral') {
                    const radii = [coil_r_max_mm];
                    for (let i = 1; i < N; i++) radii.push(radii[i - 1] - p_mm[i - 1]);
                    if(radii.some(r => r <= 0)) throw new Error("Pitch is too large, results in non-positive coil radius.");
                    coilLength_mm = radii.reduce((sum, r) => sum + 2 * Math.PI * r, 0);
                } else {
                    coilLength_mm = Math.PI * OD_mm * N;
                }

                const coilLength_m = coilLength_mm / 1000;
                const area = Math.PI * (wire_r_mm / 1000) ** 2;
                const DCR = coilLength_m / (area * sigma);

                const r0_m = wire_r_mm / 1000;
                const p_list_m = p_mm.map(v => v / 1000);

                const Hm = Array.from({ length: N }, (_, i) => {
                    const m = i + 1;
                    return this.f_H_asy(r0_m, p_list_m, m, N) + this.f_H_sym(r0_m, p_list_m, m, N);
                });

                const mu0 = 4 * Math.PI * 1e-7;
                const sd = 1 / Math.sqrt(Math.PI * f * mu0 * sigma);
                const x = 2 * r0_m / sd;
                const sumH2 = Hm.reduce((s, h) => s + h * h, 0);

                const Gp = (r0_m > sd)
                    ? (8 * sumH2 * Math.PI ** 2 * sd ** 2 * x ** 3 * (x - 1) / (((2 * x + 1) ** 2 + 2))) / N
                    : (Math.PI ** 2 * r0_m ** 2 * sumH2 * (x / 2) ** 4 / (1 + (x / 2) ** 4 / 48)) / N;

                const Rskin = (r0_m > sd)
                    ? (0.25 + r0_m / (2 * sd) + 3 * sd / (32 * r0_m)) / (Math.PI * r0_m ** 2 * sigma) * coilLength_m
                    : (1 + (r0_m / sd) ** 4 / 48) / (Math.PI * r0_m ** 2 * sigma) * coilLength_m;

                const Rprox = Gp * Rskin;
                const Rtotal = (1 + Gp) * Rskin;

                const rout_m = coil_r_max_mm / 1000;
                let L_self = (shape === 'helical')
                    ? this.f_L_self_helical(r0_m, N, rout_m, p_list_m)
                    : this.f_L_self_spiral(r0_m, N, rout_m, p_list_m);
                const L_uH = L_self * 1e6;

                // Update UI
                this.elements.coilLength.textContent = coilLength_mm.toFixed(3);
                this.elements.DCR.textContent = DCR.toFixed(6);
                this.elements.Rskin.textContent = Rskin.toFixed(6);
                this.elements.Rprox.textContent = Rprox.toFixed(6);
                this.elements.Rtotal.textContent = Rtotal.toFixed(6);
                this.elements.Gp.textContent = Gp.toFixed(3);
                this.elements.L.textContent = L_uH.toFixed(3);

            } catch (e) {
                alert("Calculation Error: " + e.message);
                console.error(e);
            }
        },
    };

    App.init();
});

