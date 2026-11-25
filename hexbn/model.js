
const WATER = {kw:0.01,h:349.81,oh:198.3};
const parseNumber = value => {
  const cleaned = String(value).trim().replace(',', '.').replace('−', '-');
  if (!/^[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(cleaned)) return null;
  const parsed = Number(cleaned);
  return Number.isFinite(parsed) ? parsed : null;
};
function constants(values) {
  if (values.some(value => parseNumber(value) === null || Math.abs(parseNumber(value)) > 30)) return null;
  return Array.from({length:3}, (_, i) => i < values.length ? 10 ** (6 - parseNumber(values[i])) : 0);
}
function species(h, [k1,k2,k3]) {
  const weights = [h ** 3, k1 * h ** 2, k1 * k2 * h, k1 * k2 * k3];
  const den = weights.reduce((sum, v) => sum + v, 0);
  const num = weights[1] + 2 * weights[2] + 3 * weights[3];
  return {den, num, charge: num / den, fractions:weights.map(weight=>weight/den)};
}
function bisect(lo,hi,below) {
  for(let i=0;i<64;i++){
    const mid=lo+(hi-lo)/2;
    if(below(mid)) lo=mid;
    else hi=mid;
  }
  return (lo+hi)/2;
}
function hydrogen(cMolar,k,water=false) {
  const c=cMolar*1e6;
  const neutral=water?Math.sqrt(WATER.kw):0;
  if(c===0)return neutral;
  if(!water&&k[1]===0)return 2*k[0]*c/(Math.sqrt(k[0]**2+4*k[0]*c)+k[0]);
  return bisect(neutral,neutral+3*c,h=>h-(water?WATER.kw/h:0) < c*species(h,k).charge);
}
function ionicConductivity(h,c,parts,water,lambdas) {
  const oh=water?WATER.kw/h:0;
  const concentrations=lambdas.map((_,i)=>c*1e6*parts.fractions[i+1]);
  const contributions=[WATER.h*h/1000,...concentrations.map((ci,i)=>lambdas[i]*ci/1000)];
  if(water)contributions.push(WATER.oh*oh/1000);
  const pure=water?(WATER.h+WATER.oh)*Math.sqrt(WATER.kw)/1000:0;
  return {oh,concentrations,contributions,pure,value:contributions.reduce((a,b)=>a+b,0)-pure};
}
function fromHydrogen(h,k,water) {
  if(!(h>0)||!Number.isFinite(h))return null;
  const parts=species(h,k);
  const oh=water?WATER.kw/h:0;
  const c=(h-oh)/parts.charge/1e6;
  return Number.isFinite(c)&&c>0?{h,c,oh,...parts}:null;
}
function fromTheory(x,k,water,lambdas) {
  const low=water?Math.sqrt(WATER.kw):0;
  const at=h=>{
    const solution=fromHydrogen(h,k,water);
    return solution?ionicConductivity(h,solution.c,solution,water,lambdas).value:0;
  };
  let high=Math.max(1,2*x*1000/WATER.h+2*low);
  while(high<1e6&&at(high)<x)high*=2;
  if(at(high)<x)return null;
  return fromHydrogen(bisect(low,high,h=>at(h)<x),k,water);
}
function fit(points) {
  if (points.length < 2) return null;
  const n = points.length;
  const xMean = points.reduce((sum, p) => sum + p.x, 0) / n;
  const yMean = points.reduce((sum, p) => sum + p.y, 0) / n;
  const sxx = points.reduce((sum, p) => sum + (p.x - xMean) ** 2, 0);
  const sxy = points.reduce((sum, p) => sum + (p.x - xMean) * (p.y - yMean), 0);
  if (!(sxx > 0) || !Number.isFinite(sxx) || !Number.isFinite(sxy)) return null;
  const a = sxy / sxx, b = yMean - a * xMean;
  const sst = points.reduce((sum, p) => sum + (p.y - yMean) ** 2, 0);
  const sse = points.reduce((sum, p) => sum + (p.y - a * p.x - b) ** 2, 0);
  return {a,b,n,xMean,yMean,sxx,sxy,sst,sse,r2:sst > 0 ? 1 - sse / sst : null};
}
const diagnosticLabels={
  ok:'Достоверная',warning:'Подозрительная',influence:'Сильно влияющая',
  repeat:'Рекомендуется переснять',unknown:'Недостаточно данных'
};
function diagnose(point,line) {
  if(!line)return {quality:'unknown',reason:'Линейный тренд недоступен'};
  const {a,b,n,xMean,yMean,sxx,sse,sst}=line;
  const predicted=a*point.x+b;
  const residual=point.y-predicted;
  const relative=point.y!==0?100*residual/Math.abs(point.y):null;
  const leverage=Math.min(1,1/n+(point.x-xMean)**2/sxx);
  const remaining=1-leverage;
  const loocv=remaining>1e-12?residual/remaining:null;
  const deletedSse=loocv===null?null:Math.max(0,sse-residual*loocv);
  const deletedSst=sst-n/(n-1)*(point.y-yMean)**2;
  const tolerance=Number.EPSILON*Math.max(sst,sse)*16;
  const deletedR2=n>2&&deletedSse!==null&&deletedSst>tolerance?1-deletedSse/deletedSst:null;
  const student=n>3&&deletedSse>tolerance?residual/Math.sqrt(deletedSse/(n-3)*remaining):null;
  const cook=n>2&&sse>tolerance&&loocv!==null?loocv**2*leverage/(2*sse/(n-2)):null;
  const concentrationError=point.back?point.back.c-point.c:null;
  const concentrationPercent=concentrationError===null?null:100*concentrationError/point.c;
  const concentrationMagnitude=Math.abs(concentrationPercent);
  const cutoff=4/n;
  const critical=[Math.abs(student)>3?'|t| > 3':'',concentrationMagnitude>20?'Ошибка C > 20%':''].filter(Boolean);
  const warnings=[Math.abs(student)>2?'|t| > 2':'',n>=4&&leverage>cutoff?'Рычаг > 4/n':'',
    concentrationMagnitude>10?'Ошибка C > 10%':''].filter(Boolean);
  let quality='unknown',reason='Часть показателей недоступна';
  if(critical.length){quality='repeat';reason=critical.join(' · ');}
  else if(cook>cutoff){quality='influence';reason='Расстояние Кука > 4/n';}
  else if(warnings.length){quality='warning';reason=warnings.join(' · ');}
  else if(student!==null&&cook!==null&&concentrationPercent!==null){quality='ok';reason='Пороги не превышены';}
  return {quality,reason,predicted,residual,relative,leverage,loocv,deletedR2,student,cook,concentrationError,concentrationPercent};
}
function invert(y,line,k,water=false,lambdas=null) {
  if(y===null||!line||!(line.a>0)||!k)return null;
  const x=(y-line.b)/line.a;
  if(!(x>0)||!Number.isFinite(x))return null;
  const result=lambdas?fromTheory(x,k,water,lambdas):fromHydrogen(x,k,water);
  return result?{...result,x}:null;
}
