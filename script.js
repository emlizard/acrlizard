// --- 핵심 계산 로직 (수정 없음) ---
Math.agm=function(a,g){var c,d=1,e=Math.sqrt,f=50;for(c=0;c<f&&1e-15<Math.abs(a-g);c++)d=a,a=(a+g)/2,g=e(d*g);return a};Math.EllipticK=function(a){return.5*Math.PI/Math.agm(1,Math.sqrt(1-a))};Math.EllipticE=function(a){for(var g,h,c=Math.sqrt(1-a),d=1,e=c,f=50,b=.25,i=1-.5*a,j=0;1e-15<Math.abs(d-e)&&j<f;)i-=b*(d-e)*(d-e),b*=2,g=.5*(d+e),h=Math.sqrt(d*e),d=g,e=h,j++;return.5*Math.PI*i/d};function f_L_self_helical(r,N,o,t){const u=4*Math.PI*1e-7,e=[0];for(let a=1;a<N;a++)e[a]=e[a-1]+t[a-1];const M=Array.from({length:N},()=>Array(N).fill(0));for(let a=0;a<N;a++)for(let l=0;l<N;l++)if(a===l)M[a][l]=u*o*(Math.log(8*o/r)-1.75);else{const i=e[a]-e[l];let s=4*o*o/(4*o*o+i*i);s=Math.min(Math.max(s,1e-12),1-1e-12);const c=Math.EllipticK(s),h=Math.EllipticE(s);M[a][l]=u*Math.sqrt(o*o/s)*((2-s)*c-2*h)}let L=0;for(let a=0;a<N;a++)for(let l=0;l<N;l++)L+=M[a][l];return L}function f_L_self_spiral(r,N,o,t){const u=4*Math.PI*1e-7,e=[o];for(let a=1;a<N;a++)e[a]=o-t.slice(0,a).reduce((s,x)=>s+x,0);const M=Array.from({length:N},()=>Array(N).fill(0));for(let a=0;a<N;a++)for(let l=0;l<N;l++)if(a===l)M[a][l]=u*e[a]*(Math.log(8*e[a]/r)-1.75);else{let s=4*e[a]*e[l]/((e[a]+e[l])*(e[a]+e[l]));s=Math.min(Math.max(s,1e-12),1-1e-12);const c=Math.EllipticK(s),h=Math.EllipticE(s);M[a][l]=u*Math.sqrt(e[a]*e[l]/s)*((2-s)*c-2*h)}let L=0;for(let a=0;a<N;a++)for(let l=0;l<N;l++)L+=M[a][l];return L}function f_H_sym(R,p,m,N){const x=N%2,NN=0===x?N/2:(N+1)/2;if(1===m||m===N)return 0;let pp_left=[],pp_right=[];if(m<=NN){for(let k=m-1,L=0;k>=1;k--)L+=p[k-1],pp_left.push(L);for(let k=m,R=0;k<m+pp_left.length;k++)R+=p[k-1],pp_right.push(R)}else{const lower=2*m-N;for(let k=m-1,L=0;k>=lower;k--)L+=p[k-1],pp_left.push(L);for(let k=m,R=0;k<m+pp_left.length;k++)R+=p[k-1],pp_right.push(R)}return pp_left.reduce((H,_,i)=>{const x1=pp_left[i]/R,x2=pp_right[i]/R,term=Math.sqrt((1+x1*x1)/Math.pow(x1*x1-1,2)+(1+x2*x2)/Math.pow(x2*x2-1,2)-2*(x1*x2-1)/((x1*x1-1)*(x2*x2-1)))/(2*Math.PI*R);return H+term},0)}function f_H_asy(R,p,m,N){const x=N%2,NN=0===x?N/2:(N+1)/2;if(m===NN&&1===x)return 0;let H=0;if(m<=NN){let offset=0;for(let k=m;k<m+m-1;k++)offset+=p[k-1];let pp_right=[];for(let k=2*m-1,R=0;k<=N-1;k++)R+=p[k-1],pp_right.push(R+offset);H=pp_right.reduce((s,v)=>s+v/(v*v+R*R)/(2*Math.PI),0)}else if(m<N){let offset=0;for(let k=m-1;k>=2*m-N;k--)offset+=p[k-1];let pp_left=[];for(let k=2*m-N-1,L=0;k>=1;k--)L+=p[k-1],pp_left.push(L+offset);H=pp_left.reduce((s,v)=>s+v/(v*v+R*R)/(2*Math.PI),0)}else{let L=0,pp_left=[];for(let k=2*m-N-1;k>=1;k--)L+=p[k-1],pp_left.push(L);H=pp_left.reduce((s,v)=>s+v/(v*v+R*R)/(2*Math.PI),0)}return H}

// --- UI 및 애플리케이션 로직 ---
document.addEventListener('DOMContentLoaded', () => {

    const themeToggleBtn = document.getElementById('theme-toggle-btn');
    const sunIcon = `<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>`;
    const moonIcon = `<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>`;

    // 테마 설정 함수
    const applyTheme = (theme) => {
        document.documentElement.setAttribute('data-theme', theme);
        themeToggleBtn.innerHTML = theme === 'light' ? moonIcon : sunIcon;
        localStorage.setItem('theme', theme);
    };

    // 테마 변경 버튼 이벤트 리스너
    themeToggleBtn.addEventListener('click', () => {
        const newTheme = document.documentElement.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
        applyTheme(newTheme);
    });
    
    // 로컬 저장소 또는 시스템 설정에 따라 초기 테마 설정
    const savedTheme = localStorage.getItem('theme') || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    applyTheme(savedTheme);

    // 페이지 로드 시 초기 계산 실행
    calculate();
});

function togglePitch() {
    const type = document.querySelector('input[name="pitchType"]:checked').value;
    document.getElementById('identicalPitchInput').style.display = type === 'identical' ? 'block' : 'none';
    document.getElementById('unequalPitchInput').style.display = type === 'unequal' ? 'block' : 'none';
}

// 임시 에러 메시지 표시 함수
function showError(message) {
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-message';
    errorDiv.textContent = message;
    document.body.appendChild(errorDiv);
    setTimeout(() => errorDiv.remove(), 4000);
}

// 메인 계산 함수
async function calculate() {
    const outputCard = document.getElementById('output-card');
    outputCard.classList.add('calculating');

    // 로딩 스피너가 UI에 표시될 시간을 벌어줌
    await new Promise(resolve => setTimeout(resolve, 50));
    
    try {
        // --- 1. 입력 값 가져오기 및 유효성 검사 ---
        const f = parseFloat(document.getElementById('f').value) * 1e6;
        const sigma = parseFloat(document.getElementById('s').value) * 1e7;
        const r0_mm = parseFloat(document.getElementById('r0').value);
        const N = parseInt(document.getElementById('N').value);
        const OD_mm = parseFloat(document.getElementById('OD').value);
        const shape = document.querySelector('input[name="coilShape"]:checked').value;
        const pitchType = document.querySelector('input[name="pitchType"]:checked').value;

        if ([f, sigma, r0_mm, N, OD_mm].some(isNaN)) {
            throw new Error("모든 숫자 입력 필드는 유효한 숫자여야 합니다.");
        }
        if (N < 1) throw new Error("코일 턴 수는 최소 1 이상이어야 합니다.");

        let p_mm;
        if (pitchType === 'identical') {
            const p0 = parseFloat(document.getElementById('p0').value);
            if (isNaN(p0)) throw new Error("피치는 유효한 숫자여야 합니다.");
            p_mm = Array(N > 1 ? N - 1 : 0).fill(p0);
        } else {
            p_mm = document.getElementById('pList').value.split(',').map(v => parseFloat(v.trim()));
            if (p_mm.some(isNaN)) throw new Error("쉼표로 구분된 모든 피치는 유효한 숫자여야 합니다.");
            if (N > 1 && p_mm.length !== N - 1) {
                 throw new Error(`${N}턴의 경우 ${N - 1}개의 피치 값을 입력해야 합니다.`);
            }
        }
        
        // --- 2. 계산 수행 ---
        const r0_m = r0_mm / 1000;
        const rout_mm = OD_mm / 2;
        const rout_m = rout_mm / 1000;
        const p_m = p_mm.map(v => v / 1000);

        let coilLength_mm;
        if (shape === 'spiral') {
            const radii = [rout_mm];
            for (let i = 1; i < N; i++) radii.push(radii[i - 1] - p_mm[i - 1]);
            if (radii.some(r => r <= 0)) throw new Error("피치가 너무 큽니다. 내부 반경이 0 또는 음수가 됩니다.");
            coilLength_mm = radii.reduce((sum, r) => sum + 2 * Math.PI * r, 0);
        } else {
            coilLength_mm = Math.PI * OD_mm * N;
        }

        const coilLength_m = coilLength_mm / 1000;
        const area = Math.PI * Math.pow(r0_m, 2);
        const DCR = coilLength_m / (area * sigma);

        const mu0 = 4 * Math.PI * 1e-7;
        const sd = 1 / Math.sqrt(Math.PI * f * mu0 * sigma);
        const x = 2 * r0_m / sd;
        
        const Hm = [];
        for (let m=1; m<=N; m++) {
            Hm.push(f_H_asy(r0_m, p_m, m, N) + f_H_sym(r0_m, p_m, m, N));
        }
        const sumH2 = Hm.reduce((s, h) => s + h * h, 0);

        const Gp = (r0_m - sd > 0)
            ? (8 * sumH2 * Math.PI**2 * sd**2 * x**3 * (x-1) / (((2*x+1)**2 + 2))) / N
            : (Math.PI**2 * r0_m**2 * sumH2 * (x/2)**4 / (1 + (x/2)**4/48)) / N;

        const Rskin = (r0_m - sd > 0)
            ? (0.25 + r0_m/(2*sd) + 3*sd/(32*r0_m)) / (Math.PI * r0_m**2 * sigma) * coilLength_m
            : (1 + Math.pow(r0_m/sd,4)/48) / (Math.PI * r0_m**2 * sigma) * coilLength_m;
        
        const Rtotal = Rskin * (1 + Gp);
        const Rprox = Rtotal - Rskin;

        let L_self = (shape === 'helical')
            ? f_L_self_helical(r0_m, N, rout_m, p_m)
            : f_L_self_spiral(r0_m, N, rout_m, p_m);
        const L_uH = L_self * 1e6;

        // --- 3. 결과 표시 ---
        // script.js 파일의 calculate 함수 내 '3. 결과 표시' 부분을 아래 코드로 교체하세요.

// --- 3. 결과 표시 ---
const results = {
    "Self-Inductance L (µH)": L_uH.toFixed(3),
    "Total AC Resistance R_total (Ω)": Rtotal.toFixed(5),
    "Wire Length ℓ (mm)": coilLength_mm.toFixed(2),
    "DC Resistance DCR (Ω)": DCR.toFixed(5),
    "Skin Effect Resistance R_skin (Ω)": Rskin.toFixed(5),
    "Proximity Effect Resistance R_prox (Ω)": Rprox.toFixed(5),
    "Proximity Factor Gp": Gp.toFixed(5)
};

const container = document.getElementById('result-container');
container.innerHTML = `<div class="result-grid"></div>`; // 이전 결과 초기화
const grid = container.querySelector('.result-grid');

// ✅ 하이라이트할 키워드 목록
const highlightKeys = ["Self-Inductance", "Total AC Resistance"];

Object.entries(results).forEach(([label, value]) => {
    const item = document.createElement('div');
    
    // ✅ 기본 클래스 설정
    item.className = 'result-item';

    // ✅ 하이라이트할 키워드가 포함된 경우, 강조 클래스 추가
    if (highlightKeys.some(key => label.includes(key))) {
        item.classList.add('result-item--highlight');
    }

    item.innerHTML = `
        <div class="result-label">${label}</div>
        <div class="result-value">${value}</div>
    `;
    grid.appendChild(item);
    
    // 업데이트 애니메이션 적용
    setTimeout(() => item.querySelector('.result-value').classList.add('updated'), 50);
});

        const container = document.getElementById('result-container');
        container.innerHTML = `<div class="result-grid"></div>`; // 이전 결과 초기화
        const grid = container.querySelector('.result-grid');
        
        Object.entries(results).forEach(([label, value]) => {
            const item = document.createElement('div');
            item.className = 'result-item';
            item.innerHTML = `
                <div class="result-label">${label}</div>
                <div class="result-value">${value}</div>
            `;
            grid.appendChild(item);
            // 업데이트 애니메이션 적용
            setTimeout(() => item.querySelector('.result-value').classList.add('updated'), 50);
        });
        
    } catch (error) {
        console.error("계산 오류:", error);
        showError(error.message);
        const container = document.getElementById('result-container');
        container.innerHTML = `<div class="placeholder-state"><p>오류: ${error.message}</p></div>`;
    } finally {
        outputCard.classList.remove('calculating');
    }
}
