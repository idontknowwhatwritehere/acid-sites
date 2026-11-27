
/* Calculation model */
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
/* View */
const seed = JSON.parse(document.getElementById('seed').textContent);
const names = Object.keys(seed);
const sources = JSON.parse(document.getElementById('sources').textContent);
const ionData=JSON.parse(document.getElementById('ion-data').textContent);
const clone = value => JSON.parse(JSON.stringify(value));
const storageKey = 'ibu-bn-calibration-v1';
const defaultBlank = 1.7;
const instrumentFactor = 1.048729;
const instrumentKey = 'ibu-bn-instrument-correction-cancelled-v1';
let instrumentCorrectionCancelled = false;
try { instrumentCorrectionCancelled = localStorage.getItem(instrumentKey)==='true'; } catch {}
const tenthsReading = new Intl.NumberFormat('en-US',{maximumFractionDigits:1,useGrouping:false});
const wholeReading = new Intl.NumberFormat('en-US',{maximumFractionDigits:0,useGrouping:false});
function instrumentReading(value,alreadyCorrected=false) {
  if(value===null||value<0)return value;
  // Remove binary floating-point noise before rounding the scaled reading.
  const measured=instrumentCorrectionCancelled||alreadyCorrected?value:Number((value*instrumentFactor).toPrecision(15));
  if(!Number.isFinite(measured))return null;
  return roundReading(measured);
}
function roundReading(value) {
  return Number((Math.abs(value)<100?tenthsReading:wholeReading).format(value));
}
function formatReading(value) {
  return Number.isFinite(value)?value.toFixed(value<100?1:0).replace('.',','):'—';
}
function displayReading(value,alreadyCorrected) {
  const raw=parseNumber(value);
  const measured=instrumentReading(raw,alreadyCorrected);
  return raw===null||raw<0||measured===null?inputText(value):formatReading(measured);
}
function inputReading(value,alreadyCorrected) {
  const reading=parseNumber(value);
  // Store readings on the original scale so both display modes recover the same measurement.
  return reading===null||instrumentCorrectionCancelled||alreadyCorrected?value:String(reading/instrumentFactor);
}
const records = Object.assign(Object.create(null),clone(seed));
let stored = {};
try { stored = JSON.parse(localStorage.getItem(storageKey) || '{}') || {}; } catch {}
for (const [name,candidate] of Object.entries(stored)) {
  if (!candidate || !Array.isArray(candidate.rows) ||
      (Object.hasOwn(seed,name) && candidate.rows.length<(name==='Адипиновая'?9:seed[name].rows.length)) ||
      !candidate.rows.every(row=>Array.isArray(row)&&row.length===2) ||
      !Array.isArray(candidate.pka) || candidate.pka.length<1 || candidate.pka.length>3) continue;
  if (!Object.hasOwn(seed,name)) names.push(name);
  if(name==='Адипиновая'&&!Array.isArray(candidate.instrumentReady)){
    const adipic=seed[name];
    candidate.instrumentReady=[];
    const updated=adipic.rows.filter(([c])=>c>=0.001);
    for(const row of updated){
      const existing=candidate.rows.findIndex(([c])=>parseNumber(c)===row[0]);
      const index=existing<0?candidate.rows.length:existing;
      candidate.rows[index]=[...row];
      const seedIndex=adipic.rows.indexOf(row);
      if(adipic.instrumentReady.includes(seedIndex))candidate.instrumentReady.push(index);
    }
  }
  records[name]=candidate;
}
let active = 'Адипиновая';
let selected = {kind:'h',row:0,step:0};
let drillOrigin = null;
let probeText = '';
let advanced = false;
let model;
const el = id => document.getElementById(id);
const fmt = (value, digits=5) => Number.isFinite(value) ? new Intl.NumberFormat('ru-RU',{maximumSignificantDigits:digits}).format(value) : '—';
const sci = value => Number.isFinite(value) ? value.toExponential(4).replace('.', ',').replace('e-', 'e−') : '—';
const inputText = value => String(value ?? '').replace('.', ',');
const make = (tag, props={}) => Object.assign(document.createElement(tag),props);
const acidTitle = name => Object.hasOwn(seed,name) ? name+' кислота' : name;
const newAcid = steps => ({pka:Array(steps).fill(''),rows:[['','']]});
const ionEntries = name => Object.hasOwn(ionData.ions,name) ? ionData.ions[name] :
  records[name].pka.map((_,i)=>['Анион '+(i+1)+' (заряд −'+(i+1)+')',null,null]);
const tabColor = name => /^#[0-9a-f]{6}$/i.test(records[name].color) ? records[name].color : '';
let colorAcid;
function persist() {
  try {
    localStorage.setItem(storageKey, JSON.stringify(records));
    return true;
  } catch {
    return false;
  }
}
function compute(name) {
  const record = records[name];
  const blankText = String(record.blank ?? defaultBlank).trim();
  const blankValue = parseNumber(blankText);
  const blankInvalid = blankText !== '' && (blankValue === null || blankValue < 0);
  const hasBlank = blankText !== '' && !blankInvalid;
  const blank = hasBlank ? blankValue : 0;
  const k = constants(record.pka);
  const water=advanced&&el('water-mode').checked;
  const wantIons=advanced&&el('ion-mode').checked;
  const entries=ionEntries(name);
  const lambdaValues=entries.map((entry,i)=>parseNumber(record.lambda?.[i] ?? (entry[1] ?? '')));
  const missing=entries.filter((_,i)=>!(lambdaValues[i]>0));
  const lambdas=wantIons&&!missing.length?lambdaValues:null;
  const rows = record.rows.map(([conc, reading], i) => {
    const c = parseNumber(conc), raw = parseNumber(reading);
    const measured=instrumentReading(raw,record.instrumentReady?.includes(i));
    const solved=k&&c!==null&&c>0?hydrogen(c,k,water):null;
    const h=Number.isFinite(solved)&&solved>0?solved:null;
    const theory=h!==null&&lambdas?ionicConductivity(h,c,species(h,k),water,lambdas):null;
    const x=theory?theory.value:h;
    return {i,c,raw,measured,y:measured===null?null:measured-blank,h,x,theory};
  });
  const points = rows.filter(row => row.x !== null && Number.isFinite(row.x) && row.measured !== null && row.measured >= 0);
  const line = fit(points);
  for (const row of rows) row.back = row.measured !== null && row.measured >= 0 ? invert(row.y,line,k,water,lambdas) : null;
  for (const point of points) point.diagnostics=diagnose(point,line);
  return {k,rows,points,line,blank,hasBlank,blankInvalid,water,wantIons,lambdas,missing,entries};
}
// Raw HexBn conductivity extrema across all source series, in µS/cm.
const hexBnRanges={"Муравьиная":[31.9,167.1],"Уксусная":[17,63.2],"Пропионовая":[21,87.3],"Яблочная":[22.9,120.4],"Щавелевая":[61.1,521],"Миндальная":[63.5,347],"Лимонная":[18.62,217],"Винная":[24.8,125],"Малоновая":[46.3,497],"Глутаровая":[32.3,114.6],"Янтарная":[55.5,132.6],"Адипиновая":[42.2,121.8]};
function appendHexBnStandards() {
  const allowed=[0.001,0.0015,0.002,0.0025,0.003,0.004,0.005,0.006,0.007,0.0075,0.008,0.009];
  for(const [name,[requiredMin,requiredMax]] of Object.entries(hexBnRanges)){
    const record=records[name];
    if(!record||record.hexBnPlanVersion===5)continue;
    const {points,line,blank,k,water,lambdas}=compute(name);
    if(!line||!(line.a>0))continue;
    if(record.hexBnPlanned?.length){
      const corrected=new Set(record.instrumentReady);
      const retained=record.rows.flatMap((row,index)=>
        record.hexBnPlanned.includes(parseNumber(row[0]))&&String(row[1]).trim()===''
          ?[]:[{row,corrected:corrected.has(index)}]);
      record.rows=retained.map(entry=>entry.row);
      if(record.instrumentReady)record.instrumentReady=retained.flatMap((entry,index)=>
        entry.corrected?[index]:[]);
    }
    record.hexBnPlanVersion=5;
    record.hexBnPlanned=[];
    const measured=points.map(point=>point.measured);
    const min=Math.min(...measured), max=Math.max(...measured);
    const occupied=new Set(record.rows.map(([c])=>parseNumber(c)));
    const expected=concentration=>{
      const h=hydrogen(concentration,k,water);
      const x=lambdas?ionicConductivity(h,concentration,species(h,k),water,lambdas).value:h;
      return line.a*x+line.b+blank;
    };
    const options=allowed.filter(c=>!occupied.has(c)).map(c=>({c,reading:expected(c)}));
    function addStandards(standards) {
      for(const {c} of standards){
        record.hexBnPlanned.push(c);
        record.rows.push([c,'']);
      }
    }
    if(requiredMin<min){
      const lower=options.filter(({reading})=>reading<min);
      const edge=lower.findLastIndex(({reading})=>reading<=requiredMin*0.95);
      if(edge>=0)addStandards(lower.slice(edge,edge+3));
    }
    if(requiredMax>max){
      const upper=options.filter(({reading})=>reading>max);
      const edge=upper.findIndex(({reading})=>reading>=requiredMax*1.05);
      if(edge>=0){
        const last=Math.min(upper.length-1,Math.max(2,edge));
        const middle=Math.ceil(last/2);
        addStandards(upper.filter((_,index)=>index===0||index===middle||index===last));
      }
    }
  }
}
const greenThreshold=0.999;
function tabShade(scores,r2,light,dark) {
  const ratio=scores.length===1?1:scores.indexOf(r2)/(scores.length-1);
  return '#'+light.map((channel,index)=>
    Math.round(channel+(dark[index]-channel)*ratio).toString(16).padStart(2,'0')).join('');
}
function hasRepeat(result) {
  return result.points.some(point=>point.diagnostics.quality==='repeat');
}
function tabStatus(result,greenScores,redScores) {
  const r2=result.line?.r2;
  if(!Number.isFinite(r2))return ['#d9e2eb','Недостаточно данных для оценки R²'];
  const score='R² '+fmt(r2,9)+' · ';
  if(r2>greenThreshold){
    const shade=tabShade(greenScores,r2,[211,245,224],[25,124,70]);
    return [shade,score+'зелёный по сравнению с другими кислотами'];
  }
  if(hasRepeat(result)){
    const shade=tabShade(redScores,r2,[255,193,200],[158,30,48]);
    return [shade,score+'Есть точки для повторного измерения'];
  }
  const statuses=new Set(result.points.map(point=>point.diagnostics.quality));
  if(statuses.has('influence'))return ['#ffca9e',score+'Есть сильно влияющие точки'];
  return ['#ffe08b',score+(statuses.has('warning')?'Есть точки для проверки':'R² не выше 0,999')];
}
function paintTabs(calibrations) {
  const greenScores=[...new Set(calibrations.map(({result})=>result.line?.r2)
    .filter(r2=>Number.isFinite(r2)&&r2>greenThreshold))].sort((a,b)=>a-b);
  const redScores=[...new Set(calibrations.filter(({result})=>hasRepeat(result))
    .map(({result})=>result.line?.r2)
    .filter(r2=>Number.isFinite(r2)&&r2<=greenThreshold))].sort((a,b)=>b-a);
  for(const {name,index,result} of calibrations){
    paintTab(el('tabs').children[index],name,tabStatus(result,greenScores,redScores));
  }
}
function paintTab(item,name,[automatic,status]) {
  const color=tabColor(name)||automatic;
  const rgb=color.slice(1).match(/../g).map(hex=>{
    const channel=parseInt(hex,16)/255;
    return channel<=0.04045?channel/12.92:((channel+0.055)/1.055)**2.4;
  });
  const luminance=0.2126*rgb[0]+0.7152*rgb[1]+0.0722*rgb[2];
  item.classList.add('custom-color');
  item.style.setProperty('--tab-color',color);
  item.dataset.autoColor=automatic;
  item.style.setProperty('--tab-ink',luminance>0.179?'#000':'#fff');
  const tab=item.querySelector('.tab');
  tab.title=status;
  tab.setAttribute('aria-description',status);
}
function showTabs() {
  const scrollLeft=el('tabs').scrollLeft;
  el('tabs').replaceChildren(...names.map(name => {
    const item=make('div',{className:'acid-tab'});
    item.setAttribute('role','presentation');
    const node = make('button',{className:'tab',type:'button',textContent:name});
    node.setAttribute('role','tab');
    node.setAttribute('aria-selected',String(name === active));
    node.onclick = () => { active=name; selected={kind:'h',row:0,step:0}; drillOrigin=null; mount(); focusAcid(name); };
    const menu=make('button',{className:'tab-menu',type:'button',textContent:'⋯',title:'Настройки вкладки'});
    menu.setAttribute('aria-label','Настройки вкладки «'+name+'»');
    menu.setAttribute('aria-haspopup','dialog');
    menu.onclick=()=>{
      colorAcid=name;
      el('color-title').textContent='Настройки вкладки «'+name+'»';
      el('tab-color').value=tabColor(name)||item.dataset.autoColor;
      el('reset-color').disabled=!tabColor(name);
      el('delete-acid').hidden=Object.hasOwn(seed,name);
      el('color-dialog').showModal();
    };
    item.append(node,menu);
    return item;
  }));
  el('tabs').scrollLeft=scrollLeft;
}
function focusAcid(name, selector='.tab') {
  el('tabs').children[names.indexOf(name)].querySelector(selector).focus();
}
function showInputs() {
  const record = records[active];
  el('sheet-title').textContent = acidTitle(active);
  const source = Object.hasOwn(sources,active) ? sources[active] : null;
  el('pk-source').hidden=!source;
  if(source) {
    el('pk-source').href=source.url;
    el('pk-source').textContent='Источник pKₐ · '+source.label;
  }
  el('range-source').textContent=source?'Лист и ячейки · '+source.range:'Пользовательские данные';
  const entries=ionEntries(active);
  el('ion-fields').replaceChildren(...entries.map(([name,value,ref,note],i)=>{
    const label=make('label',{textContent:name});
    const input=make('input',{id:'lambda-'+i,type:'text',inputMode:'decimal',value:inputText(record.lambda?.[i]??value??''),placeholder:'нет λ°'});
    input.setAttribute('aria-label','λ° '+name+', См·см²/моль');
    if(note)input.title=note;
    input.oninput=()=>{
      record.lambda={...record.lambda,[i]:input.value};
      persist();refresh();
    };
    label.append(input);
    return label;
  }));
  el('ion-refs').replaceChildren('См·см²/моль · ',...Array.from(new Set(entries.map(e=>e[2]).filter(Boolean))).flatMap(ref=>{
    const [label,url]=ionData.refs[ref];
    return [make('a',{href:url,target:'_blank',rel:'noopener noreferrer',textContent:label}),' · '];
  }),'H⁺ и OH⁻: NIST');
  el('blank').value = inputText(record.blank ?? defaultBlank);
  el('blank').onfocus = () => select('blank',selected.row);
  el('blank').oninput = () => {record.blank=el('blank').value; persist(); refresh();};
  el('constants').replaceChildren(...record.pka.map((value, step) => {
    const wrap = make('span',{className:'constant'});
    const label = make('label',{htmlFor:'pka-'+step,textContent:'pKₐ'+'₁₂₃'[step]});
    const input = make('input',{id:'pka-'+step,value:inputText(value),inputMode:'decimal'});
    input.dataset.pka = step;
    input.onfocus = () => select('ka',selected.row,step);
    input.oninput = () => { record.pka[step]=input.value; persist(); refresh(); };
    const ka = make('button',{id:'ka-'+step,type:'button'});
    ka.onclick = () => select('ka',selected.row,step);
    wrap.append(label,input,ka);
    return wrap;
  }));
  el('save-status').textContent='';
  el('rows').replaceChildren(...record.rows.map((row, i) => {
    const alreadyCorrected=record.instrumentReady?.includes(i);
    const tr = make('tr');
    const number = make('th',{scope:'row'});
    const numberContent=make('span',{className:'row-number',textContent:i+1});
    if(!Object.hasOwn(seed,active)||i>=seed[active].rows.length){
      const remove=make('button',{type:'button',className:'remove-row',textContent:'×',title:'Удалить строку'});
      remove.setAttribute('aria-label','Удалить строку '+(i+1));
      remove.onclick=()=>{
        if(row.some(value=>String(value).trim())&&!confirm('Удалить строку '+(i+1)+' и её измерение?'))return;
        record.rows.splice(i,1);
        if(record.instrumentReady)record.instrumentReady=record.instrumentReady.filter(index=>index!==i)
          .map(index=>index>i?index-1:index);
        persist();
        drillOrigin=null;
        const next=Math.max(0,Math.min(i,record.rows.length-1));
        selected={kind:'h',row:next,step:0};
        showInputs();
        refresh();
        (el('conc-'+next)||el('add-row')).focus();
      };
      numberContent.append(remove);
    }
    number.append(numberContent);
    tr.append(number);
    ['conc','reading','h','back','theory'].forEach((kind, col) => {
      const td = make('td');
      const control = col < 2
        ? make('input',{value:col===1?displayReading(row[col],alreadyCorrected):inputText(row[col]),type:'text',inputMode:'decimal',spellcheck:false})
        : make('button',{type:'button',className:'value'});
      control.id = kind+'-'+i;
      control.setAttribute('aria-label','Стандарт '+(i+1)+' · '+['концентрация','показание','расчёт H⁺','концентрация по тренду','теоретическая проводимость'][col]);
      control.onfocus = () => select(kind,i);
      if (col < 2) control.oninput = () => { row[col]=col===1?inputReading(control.value,alreadyCorrected):control.value; persist(); refresh(); };
      else control.onclick = () => select(kind,i);
      if(kind==='reading')control.onblur=()=>{control.value=displayReading(row[col],alreadyCorrected);};
      td.append(control);
      tr.append(td);
    });
    return tr;
  }));
}
function select(kind,row=selected.row,step=selected.step) {
  selected={kind,row,step};
  paintSelection();
  showFormula();
  paintPoint();
}
function paintSelection() {
  document.querySelectorAll('.selected').forEach(node => node.classList.remove('selected'));
  const id = selected.kind==='ka' ? 'ka-'+selected.step : selected.kind==='blank' ? 'blank' : selected.kind+'-'+selected.row;
  el(id)?.classList.add('selected');
}
function paintPoint() {
  el('chart').querySelectorAll('[data-point]').forEach(node => node.classList.toggle('active',
    selected.kind !== 'probe' && Number(node.dataset.point) === selected.row));
}
function updateMainFormula({water,lambdas,wantIons,missing}) {
  const useIons=Boolean(lambdas);
  document.querySelectorAll('.inverse-h').forEach(term=>term.toggleAttribute('hidden',water));
  document.querySelectorAll('.inverse-h-water').forEach(term=>term.toggleAttribute('hidden',!water));
  el('main-conductivity').hidden=!useIons;
  el('main-hydroxide-term').toggleAttribute('hidden',!water||!useIons);
  el('main-theory-baseline').toggleAttribute('hidden',!water||!useIons);
  el('main-mode').textContent=useIons?'Ось X: расчётная κ, мкСм/см.':
    wantIons&&missing.length?'Ось X: [H⁺]. Для режима Кольрауша нужны λ°: '+missing.map(entry=>entry[0]).join(', ')+'.':
    'Ось X: [H⁺], мкмоль/л.';
  el('main-conductivity-note').textContent='cᵢ — концентрация каждой анионной формы кислоты, λᵢ — её предельная молярная проводимость. Концентрации подставляются в мкмоль/л.'+
    (water?' Вклад чистой воды в модели равен 0,054811 мкСм/см.':'');
}
function refresh() {
  model=compute(active);
  try {
    localStorage.setItem('ibu-bn-model-v1',JSON.stringify({water:model.water,ions:model.wantIons}));
  } catch {}
  const calibrations=names.map((name,index)=>({name,index,result:name===active?model:compute(name)}));
  paintTabs(calibrations);
  updateMainFormula(model);
  records[active].pka.forEach((value, i) => {
    const parsed=parseNumber(value);
    el('pka-'+i).setAttribute('aria-invalid',String(parsed===null || Math.abs(parsed)>30));
    el('ka-'+i).textContent='Kₐ '+fmt(model.k?.[i],4);
  });
  for(const row of model.rows) {
    const conc=el('conc-'+row.i), reading=el('reading-'+row.i);
    conc.setAttribute('aria-invalid',String(conc.value.trim()!=='' && !(row.c>0)));
    reading.setAttribute('aria-invalid',String(reading.value.trim()!=='' && (row.raw===null || row.raw<0)));
    const planned=records[active].hexBnPlanned?.includes(row.c)&&row.raw===null;
    const expected=planned&&model.line&&row.x!==null?model.line.a*row.x+model.line.b+model.blank:null;
    reading.placeholder=Number.isFinite(expected)?'Ожид. '+formatReading(expected):'';
    reading.title=planned?'План HexBn. Введите измеренную κ, мкСм/см.':'';
    el('h-'+row.i).textContent=fmt(row.h,6);
    el('back-'+row.i).textContent=sci(row.back?.c);
    el('theory-'+row.i).textContent=fmt(row.theory?.value,6);
    el('theory-'+row.i).parentElement.hidden=!model.lambdas;
  }
  for(const id of ['theory-col','theory-letter','theory-head','theory-option']) el(id).hidden=!model.lambdas;
  el('balance-option').hidden=!model.water;
  if ((selected.kind==='theory'&&!model.lambdas)||(selected.kind==='balance'&&!model.water)) {
    selected={...selected,kind:'h'};
    paintSelection();
  }
  if(model.wantIons&&model.missing.length)el('ion-editor').open=true;
  el('model-status').textContent=model.wantIons&&model.missing.length?
    'Ось X: [H⁺]. Для суммы ионов укажи λ°: '+model.missing.map(e=>e[0]).join(', ')+'.':
    model.lambdas?'Ось X: расчётная проводимость всех ионов.':'Ось X: расчётная [H⁺].';
  el('blank').setAttribute('aria-invalid',String(model.blankInvalid));
  el('blank-state').textContent=model.blankInvalid?'Введи неотрицательное число.':
    model.hasBlank?'Вычитается из стандартов и пробы.':'Фон не задан.';
  el('notice').textContent=!model.k?'Проверь pKa: для каждой ступени нужно число.':
    model.points.length < model.rows.length ? 'Строки с пустыми или некорректными данными исключены из тренда.' :
    model.line && !(model.line.a>0) ? 'Обратный пересчёт требует положительного наклона тренда.' : '';
  renderChart();
  renderAdditionalCharts(calibrations);
  showFormula();
}
function mount() {
  showTabs(); showInputs(); refresh(); paintSelection();
}
function axisCeiling(value) {
  if (!(value > 0)) return 1;
  const base=10 ** Math.floor(Math.log10(value/5));
  const step=[1,2,2.5,5,10].map(m=>m*base).find(v=>v>=value/5) || 10*base;
  return Math.ceil(value/step) * step;
}
function hidePointTooltip() {
  const tooltip=el('point-tooltip');
  tooltip.hidden=true;
  document.querySelector('[aria-describedby="point-tooltip"]')?.removeAttribute('aria-describedby');
}
function showPointTooltip(node,point) {
  const {line,points,hasBlank,lambdas}=model;
  const d=point.diagnostics;
  const units=' мкСм/см';
  const quantity=(value,unit='')=>Number.isFinite(value)?fmt(value,6)+unit:'—';
  const values=[
    ['C стандарта',sci(point.c)+' моль/л'],
    ['[H⁺]',fmt(point.h,6)+' мкмоль/л'],
    ...(lambdas ? [['κ расчётная (X)',fmt(point.x,6)+units]] : []),
    ['Показание κ',formatReading(point.measured)+units],
    ...(hasBlank ? [['κ − κ₀ (Y)',fmt(point.y,6)+units]] : []),
    ...(line ? [
      ['Расчётная κ по тренду',quantity(d.predicted,units)],
      ['Остаток Δκ',quantity(d.residual,units)],
      ['Абсолютный остаток |Δκ|',quantity(Math.abs(d.residual),units)],
      ['Относительный остаток',quantity(d.relative,'%')],
      ['Стьюдентизированный остаток',quantity(d.student)],
      ['Расстояние Кука',quantity(d.cook)],
      ['Рычаг влияния (leverage)',quantity(d.leverage)],
      ['R² тренда',fmt(line.r2,6)],
      ['R² без точки',quantity(d.deletedR2)],
      ['Ошибка прогноза LOOCV',quantity(d.loocv,units)],
      ['RMSE тренда',fmt(Math.sqrt(line.sse/line.n),5)+units],
      ['Ошибка обратного расчёта C',Number.isFinite(d.concentrationError)?
        sci(d.concentrationError)+' моль/л · '+quantity(d.concentrationPercent,'%'):'—']
    ] : []),
    ['Точек в калибровке',String(points.length)]
  ];
  showChartTooltip(node,'Стандарт '+(point.i+1)+' · '+diagnosticLabels[d.quality],values,d.reason);
}
function showChartTooltip(node,title,values,note) {
  const tooltip=el('point-tooltip');
  hidePointTooltip();
  const details=make('dl');
  details.append(...values.flatMap(([label,value])=>[
    make('dt',{textContent:label}),make('dd',{textContent:value})
  ]));
  tooltip.replaceChildren(make('strong',{textContent:title}));
  if(note)tooltip.append(make('p',{className:'detail',textContent:note}));
  tooltip.append(details);
  tooltip.hidden=false;
  node.setAttribute('aria-describedby','point-tooltip');
  const anchor=node.getBoundingClientRect();
  const x=anchor.left+anchor.width/2;
  const left=Math.max(8,Math.min(x-tooltip.offsetWidth/2,document.documentElement.clientWidth-tooltip.offsetWidth-8));
  const above=anchor.top-tooltip.offsetHeight-8;
  const top=above>=8?above:Math.max(8,Math.min(anchor.bottom+8,innerHeight-tooltip.offsetHeight-8));
  tooltip.style.left=left+'px';
  tooltip.style.top=top+'px';
}
function renderChart() {
  const {points,line}=model;
  hidePointTooltip();
  el('metric').textContent='R² '+fmt(line?.r2,6);
  el('trend').replaceChildren(...(line ? [
    model.hasBlank?'κ − κ₀ = ':'κ = ',
    make('span',{id:'slope-value',textContent:fmt(line.a,5)}), model.lambdas?' × κтеор ':' × [H⁺] ',
    make('span',{id:'intercept-value',textContent:(line.b>=0 ? '+ ' : '− ')+fmt(Math.abs(line.b),4)})
  ] : ['Тренд недоступен']));
  if (!points.length) {
    el('chart').replaceChildren(make('p',{className:'plot-empty',textContent:'Введи концентрацию и показание для отображения точек.'}));
    return;
  }
  const width=630, height=324, left=67, right=606, top=18, bottom=272;
  const minX=Math.min(...points.map(p=>p.x)),maxX=Math.max(...points.map(p=>p.x));
  const ymax=axisCeiling(1.06*Math.max(...points.map(p=>p.y),line ? line.a*minX+line.b : 0,line ? line.a*maxX+line.b : 0));
  const ymin=Math.min(0,...points.map(p=>p.y),line ? line.a*minX+line.b : 0,line ? line.a*maxX+line.b : 0);
  const xmax=axisCeiling(maxX*1.06);
  const px=x=>left+(right-left)*x/xmax;
  const py=y=>bottom-(bottom-top)*(y-ymin)/(ymax-ymin);
  const elements=[];
  for(let i=0;i<=5;i++){
    const x=xmax*i/5,y=ymin+(ymax-ymin)*i/5;
    elements.push('<path class="grid" d="M'+px(x)+' '+top+'V'+bottom+' M'+left+' '+py(y)+'H'+right+'"/>');
    elements.push('<text x="'+px(x)+'" y="293" text-anchor="middle">'+fmt(x,4)+'</text>');
    elements.push('<text x="59" y="'+(py(y)+4)+'" text-anchor="end">'+fmt(y,4)+'</text>');
  }
  elements.push('<path class="axis" fill="none" d="M'+left+' '+top+'V'+bottom+'H'+right+'"/>');
  if(line)elements.push('<line class="trend" x1="'+px(minX)+'" y1="'+py(line.a*minX+line.b)+'" x2="'+px(maxX)+'" y2="'+py(line.a*maxX+line.b)+'"/>');
  for(const p of points) {
    const quality=p.diagnostics.quality;
    const label='Стандарт '+(p.i+1)+(model.lambdas?' · κтеор ':' · [H⁺] ')+fmt(p.x,6)+(model.hasBlank?' · κ − κ₀ ':' · κ ')+fmt(p.y,6)+' мкСм/см · '+diagnosticLabels[quality];
    elements.push('<circle class="point" tabindex="0" role="button" aria-label="'+label+'" data-point="'+p.i+'" data-quality="'+quality+'" cx="'+px(p.x)+'" cy="'+py(p.y)+'" r="5"/>');
  }
  elements.push('<text class="axis-label" x="335" y="318" text-anchor="middle">'+(model.lambdas?'Расчётная κ, мкСм/см':'[H⁺] расчётная, мкмоль/л')+'</text>');
  elements.push('<text class="axis-label" x="18" y="151" text-anchor="middle" transform="rotate(-90 18 151)">'+(model.hasBlank?'κ − κ₀':'κ')+', мкСм/см</text>');
  el('chart').innerHTML='<svg class="plot" viewBox="0 0 '+width+' '+height+'" aria-label="'+(model.lambdas?'Измеренная проводимость от расчётной':'Показание прибора от концентрации H⁺')+'">'+elements.join('')+'</svg>';
  el('chart').querySelectorAll('[data-point]').forEach(node => {
    const point=model.rows[Number(node.dataset.point)];
    const describe=()=>showPointTooltip(node,point);
    const choose=()=>{select(model.lambdas?'theory':'h',point.i); describe();};
    node.onpointerenter=describe;
    node.onpointerleave=event=>{if(event.pointerType!=='touch')hidePointTooltip();};
    node.onfocus=describe;
    node.onblur=hidePointTooltip;
    node.onclick=choose;
    node.onkeydown=event=>{
      if(event.key==='Escape')hidePointTooltip();
      if(event.key==='Enter'||event.key===' '){event.preventDefault();choose();}
    };
  });
  paintPoint();
}
const chartColors=['#1c5395','#087c71','#9f3c74','#aa530b','#6844aa','#956b00',
  '#176e97','#ac343b','#53732a','#5343a3','#7b4f34','#246b58'];

function drawAdditionalChart(id, {xLabel,yLabel,dots=[],curves=[],reference=[],logX=false,
    xRange,yRange,centerY=false,empty='Недостаточно данных для графика.'}) {
  const host=el(id+'-chart');
  const legend=el(id+'-legend');
  if(legend){
    legend.replaceChildren(...curves.map(series=>{
      const item=make('li',{textContent:series.label});
      item.style.setProperty('--series-color',series.color);
      return item;
    }));
    legend.hidden=!curves.length;
  }
  if(!dots.length&&!curves.length){
    host.replaceChildren(make('p',{className:'plot-empty',textContent:empty}));
    return;
  }
  const width=570,height=295,left=85,right=536,top=16,bottom=233;
  const all=[...dots,...curves.flatMap(series=>series.points),...reference];
  const xPosition=x=>logX?Math.log10(x):x;
  const xValues=all.map(({x})=>xPosition(x));
  const yValues=all.map(({y})=>y);
  const [xmin,xmax]=xRange ?? (logX?
    [Math.floor(Math.min(...xValues)),Math.max(Math.floor(Math.min(...xValues))+1,Math.ceil(Math.max(...xValues)))]:
    [0,axisCeiling(Math.max(...xValues)*1.03)]);
  const yPeak=axisCeiling(Math.max(...yValues.map(Math.abs))*1.06);
  const [ymin,ymax]=yRange ?? (centerY?[-yPeak,yPeak]:[Math.min(0,...yValues),yPeak]);
  const px=x=>left+(right-left)*(xPosition(x)-xmin)/(xmax-xmin);
  const py=y=>bottom-(bottom-top)*(y-ymin)/(ymax-ymin);
  const yTicks=Array.from({length:6},(_,i)=>ymin+(ymax-ymin)*i/5);
  const logStep=Math.max(1,Math.ceil((xmax-xmin)/5));
  const xTicks=logX?Array.from({length:Math.floor((xmax-xmin)/logStep)+1},(_,i)=>10**(xmin+i*logStep)):
    Array.from({length:6},(_,i)=>xmin+(xmax-xmin)*i/5);
  const svg=document.createElementNS('http://www.w3.org/2000/svg','svg');
  svg.setAttribute('class','plot');
  svg.setAttribute('viewBox','0 0 '+width+' '+height);
  svg.setAttribute('role',dots.length?'group':'img');
  svg.setAttribute('aria-labelledby',id+'-title');
  svg.setAttribute('aria-label',xLabel+'; '+yLabel);
  function shape(tag,attributes,label,parent=svg) {
    const node=document.createElementNS(svg.namespaceURI,tag);
    for(const [name,value] of Object.entries(attributes))node.setAttribute(name,value);
    if(label!==undefined)node.textContent=label;
    parent.append(node);
    return node;
  }
  const formatTick=value=>value!==0&&(Math.abs(value)<0.001||Math.abs(value)>=1e6)?
    value.toExponential(1).replace('.0e','e').replace('.',','):fmt(value,4);
  for(const x of xTicks){
    shape('path',{class:'grid',d:'M'+px(x)+' '+top+'V'+bottom});
    shape('text',{x:px(x),y:bottom+20,'text-anchor':'middle'},formatTick(x));
  }
  for(const y of yTicks){
    shape('path',{class:'grid',d:'M'+left+' '+py(y)+'H'+right});
    shape('text',{x:left-8,y:py(y)+4,'text-anchor':'end'},formatTick(y));
  }
  shape('path',{class:'axis',fill:'none',d:'M'+left+' '+top+'V'+bottom+'H'+right});
  function drawLine(points,attributes,label) {
    const path=points.map(({x,y},i)=>(i?'L':'M')+px(x)+' '+py(y)).join(' ');
    const layer=shape('path',{...attributes,d:path});
    if(label)shape('title',{},label,layer);
  }
  if(centerY)shape('line',{class:'chart-reference',x1:left,x2:right,y1:py(0),y2:py(0)});
  if(reference.length)drawLine(reference,{class:'chart-reference'});
  for(const series of curves){
    drawLine(series.points,{class:'chart-curve',stroke:series.color,
      'stroke-dasharray':series.dash ?? 'none'},series.label);
  }
  for(const dot of dots){
    const title=dot.label ?? 'Стандарт '+(dot.i+1);
    const point=shape('circle',{class:'chart-dot',fill:dot.color ?? chartColors[0],cx:px(dot.x),cy:py(dot.y),
      r:dot.radius ?? 4.5,tabindex:0,role:'button','aria-label':title+' · '+formatTick(dot.x)+'; '+fmt(dot.y,6)});
    if(dot.comparison)point.setAttribute('data-comparison','');
    const show=()=>{
      if(point.hasAttribute('aria-describedby')){
        hidePointTooltip();
        return;
      }
      const {values=[],note}=dot.describe?.() ?? standardPointDetails(model.rows[dot.i]);
      showChartTooltip(point,title,[[xLabel,fmt(dot.x,7)],[yLabel,fmt(dot.y,7)],...values],note);
    };
    point.onclick=show;
    point.onkeydown=event=>{
      if(event.key==='Enter'||event.key===' '){
        event.preventDefault();
        show();
      }
    };
  }
  shape('text',{class:'axis-label',x:(left+right)/2,y:height-11,'text-anchor':'middle'},xLabel);
  const yMiddle=(top+bottom)/2;
  shape('text',{class:'axis-label',x:19,y:yMiddle,'text-anchor':'middle',
    transform:'rotate(-90 19 '+yMiddle+')'},yLabel);
  host.replaceChildren(svg);
}

function standardPointDetails(row, diagnostics=row.diagnostics, line=model.line, back=row.back) {
  const quantity=(value,unit='')=>Number.isFinite(value)?fmt(value,6)+unit:'—';
  const values=[
    ['C стандарта',quantity(row.c,' моль/л')],
    ['Измеренная κ',Number.isFinite(row.measured)?formatReading(row.measured)+' мкСм/см':'—'],
    ['C по тренду',quantity(back?.c,' моль/л')],
  ];
  if(line&&diagnostics){
    values.push(
      ['κ по тренду',quantity(diagnostics.predicted,' мкСм/см')],
      ['Остаток Δκ',quantity(diagnostics.residual,' мкСм/см')],
      ['Относительный остаток κ',quantity(diagnostics.relative,'%')],
      ['Ошибка пересчёта ΔC',quantity(diagnostics.concentrationError,' моль/л')],
      ['Относительная ошибка C',quantity(diagnostics.concentrationPercent,'%')],
      ['Стьюдентизированный остаток',quantity(diagnostics.student)],
      ['Расстояние Кука',quantity(diagnostics.cook)],
      ['R² тренда',quantity(line.r2)],
    );
  }
  return {values,note:diagnostics?
    diagnosticLabels[diagnostics.quality]+' · '+diagnostics.reason:undefined};
}

function renderAdditionalCharts(calibrations) {
  const {k,rows,points,line,lambdas}=model;
  const concentration='C стандарта, моль/л';
  const signal='κ, мкСм/см';
  const hydrogenLabel='[H⁺] расчётная, мкмоль/л';
  const theoryLabel='κ расчётная, мкСм/см';
  const standards=rows.filter(row=>row.c>0);
  const recovered=standards.filter(row=>row.back?.c>0);
  const dot=(row,x,y)=>({i:row.i,x,y});
  drawAdditionalChart('raw',{xLabel:concentration,yLabel:signal,
    dots:standards.filter(row=>row.measured!==null&&row.measured>=0).map(row=>dot(row,row.c,row.measured))});
  drawAdditionalChart('error',{xLabel:concentration+' (лог. шкала)',yLabel:'Ошибка, %',logX:true,centerY:true,
    dots:recovered.map(row=>dot(row,row.c,100*(row.back.c/row.c-1))).filter(point=>Number.isFinite(point.y))});
  drawAdditionalChart('residual',{xLabel:lambdas?theoryLabel:hydrogenLabel,yLabel:'Остаток κ, мкСм/см',centerY:true,
    dots:line?points.map(row=>dot(row,row.x,row.y-(line.a*row.x+line.b))):[]});

  const pka=records[active].pka.map(parseNumber);
  const pHMin=k?Math.min(0,Math.floor(Math.min(...pka))-2):0;
  const pHMax=k?Math.max(14,Math.ceil(Math.max(...pka))+2):14;
  const formulas=[['HA','A⁻'],['H₂A','HA⁻','A²⁻'],['H₃A','H₂A⁻','HA²⁻','A³⁻']];
  const speciesCurves=k?formulas[pka.length-1].map((label,index)=>({
    label,color:chartColors[index],
    points:Array.from({length:141},(_,i)=>{
      const x=pHMin+(pHMax-pHMin)*i/140;
      return {x,y:species(10**(6-x),k).fractions[index]*100};
    }),
  })):[];
  const speciesLabels=formulas[pka.length-1];
  const speciesDots=k?pka.map((x,i)=>{
    const fractions=species(10**(6-x),k).fractions;
    return {x,y:fractions[i]*100,label:'pKₐ'+['₁','₂','₃'][i]+' · '+speciesLabels[i],color:chartColors[i],
      describe:()=>({values:speciesLabels.map((label,index)=>[label,fmt(fractions[index]*100,5)+'%'])})};
  }):[];
  drawAdditionalChart('species',{xLabel:'pH',yLabel:'Доля формы, %',
    xRange:[pHMin,pHMax],yRange:[0,100],curves:speciesCurves,dots:speciesDots,empty:'Введи pKₐ для каждой ступени диссоциации.'});

  const concentrations=standards.map(row=>row.c);
  const min=Math.min(...concentrations),max=Math.max(...concentrations);
  const protonsAt=x=>hydrogen(x,k,false)/(x*1e6);
  const protonCurve=k&&max>min?Array.from({length:101},(_,i)=>{
    const x=10**(Math.log10(min)+(Math.log10(max)-Math.log10(min))*i/100);
    return {x,y:protonsAt(x)};
  }):[];
  drawAdditionalChart('protons',{xLabel:'C кислоты, моль/л (лог. шкала)',yLabel:'h/C, протонов на молекулу',
    logX:true,yRange:[0,pka.length],curves:protonCurve.length?[{color:chartColors[0],points:protonCurve}]:[],
    dots:k?standards.map(row=>dot(row,row.c,protonsAt(row.c))).filter(point=>Number.isFinite(point.y)):[]});

  const recoveredMax=axisCeiling(Math.max(0,...recovered.flatMap(row=>[row.c,row.back.c]))*1.03);
  drawAdditionalChart('recovered',{xLabel:concentration,yLabel:'C по тренду, моль/л',
    xRange:[0,recoveredMax],yRange:[0,recoveredMax],reference:[{x:0,y:0},{x:recoveredMax,y:recoveredMax}],
    dots:recovered.map(row=>dot(row,row.c,row.back.c))});

  const fitted=calibrations.filter(({result})=>result.line);
  const commonTheory=fitted.length>0&&fitted.every(({result})=>result.lambdas);
  const comparison=fitted.flatMap(({name,index,result})=>{
    const refit=!commonTheory&&result.lambdas;
    const calibrationPoints=refit?result.points.map(row=>({...row,x:row.h})):result.points;
    const calibrationLine=refit?fit(calibrationPoints):result.line;
    if(!calibrationLine)return [];
    const xValues=calibrationPoints.map(row=>row.x);
    const color=chartColors[index%chartColors.length];
    return [{
      color,dash:[undefined,'7 4','2 3'][Math.floor(index/4)%3],
      label:acidTitle(name),
      points:[Math.min(...xValues),Math.max(...xValues)].map(x=>({x,y:calibrationLine.a*x+calibrationLine.b})),
      dots:calibrationPoints.map(row=>({
        x:row.x,y:row.y,color,radius:3.2,comparison:true,
        label:acidTitle(name)+' · стандарт '+(row.i+1),
        describe:()=>{
          const back=invert(row.y,calibrationLine,result.k,result.water,commonTheory?result.lambdas:null);
          return standardPointDetails(row,diagnose({...row,back},calibrationLine),calibrationLine,back);
        },
      })),
    }];
  });
  const note=el('comparison-axis-note');
  note.hidden=!model.wantIons||commonTheory||!comparison.length;
  note.textContent='Для всех кислот использована ось [H⁺], поскольку для части калибровок не заданы все λ°.';
  drawAdditionalChart('comparison',{xLabel:commonTheory?theoryLabel:hydrogenLabel,
    yLabel:'κ с учётом фона, мкСм/см',curves:comparison,dots:comparison.flatMap(series=>series.dots)});
}

function drill(kind) {
  drillOrigin = {...selected};
  select(kind);
  el('formula-back').focus({preventScroll:true});
  el('formula-title').scrollIntoView({block:'nearest',behavior:motion()?'smooth':'auto'});
}
function chargeContext(data, line, k, origin) {
  const isProbe = origin === 'probe';
  if(origin === 'h'||origin === 'balance') return data.h!==null&&k ? {h:data.h,parts:species(data.h,k),ids:['h-'+data.i]} : null;
  const measured=isProbe?instrumentReading(parseNumber(probeText)):data.measured;
  const result=measured!==null&&measured>=0?invert(measured-model.blank,line,k,model.water,model.lambdas):null;
  return result ? {h:result.h,parts:result,ids:[isProbe?'probe':'reading-'+data.i,'slope-value','intercept-value',...(model.hasBlank?['blank']:[])]} : null;
}
function showFormula() {
  const {kind,row,step}=selected;
  const index=row+1;
  const titles={ka:'Kₐ'+'₁₂₃'[step]+' из pKₐ',blank:'B'+index+' · Поправка на фон',
    reading:'B'+index+' · Поправка на фон',conc:'A'+index+' · Исходная концентрация',
    h:'C'+index+' · Расчёт [H⁺]',fit:'Линейная регрессия',r2:'Коэффициент R²',
    probe:'Концентрация пробы',back:'D'+index+' · Концентрация по тренду',
    num:'Расчёт N · числитель',den:'Расчёт D · знаменатель',theory:'E'+index+' · Проводимость по Кольраушу',balance:'Баланс зарядов и вода'};
  el('formula-back').hidden=kind!=='num'&&kind!=='den';
  el('probe').hidden=kind!=='probe';
  el('mode').value=['h','back','fit','r2','probe','blank','theory','balance'].includes(kind)?kind:'';
  el('selection').textContent=titles[kind];
  showNumeric(kind,model.rows[row]??{i:0,h:null,c:null,raw:null},model.line,model.k,step);
}
function highlightSources(ids, enabled) {
  if(enabled) document.querySelectorAll('.source-highlight').forEach(node=>node.classList.remove('source-highlight'));
  ids.forEach(id=>el(id)?.classList.toggle('source-highlight',enabled));
}
const motion=()=>!matchMedia('(prefers-reduced-motion: reduce)').matches;
function pulseSources(ids) {
  highlightSources(ids,true);
  const source=ids.map(el).find(Boolean);
  if(!source) return;
  source.scrollIntoView({block:'nearest',behavior:motion()?'smooth':'auto'});
  if(motion()) source.animate([{transform:'scale(1)'},{transform:'scale(1.075)'},{transform:'scale(1)'}],
    {duration:350,easing:'ease-out'});
}
function numberToken(value, tint, ids=[], tooltip='', scientific=false, destination=null) {
  const shown=typeof value==='string'?value:scientific?sci(value):fmt(value,7);
  const span=make('span',{className:'num-token',textContent:shown});
  span.dataset.tint=tint;
  if(ids.length){
    span.tabIndex=0;
    span.dataset.action='source';
    span.setAttribute('role','button');
    if(destination) span.dataset.drill=typeof destination==='string'?destination:'factor';
    span.title=tooltip+(destination?' · открыть расчёт':' · показать источник');
    const activate=()=>typeof destination==='function'?destination():destination?drill(destination):pulseSources(ids);
    span.addEventListener('click',activate);
    span.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();activate();}});
    span.setAttribute('aria-label',span.textContent+' · '+tooltip);
    span.addEventListener('pointerenter',()=>highlightSources(ids,true));
    span.addEventListener('pointerleave',()=>{if(document.activeElement!==span)highlightSources(ids,false);});
    span.addEventListener('focus',()=>highlightSources(ids,true));
    span.addEventListener('blur',()=>{if(!span.matches(':hover'))highlightSources(ids,false);});
  }
  if(!tooltip)return span;
  const label=make('small',{className:'number-label',textContent:tooltip});
  if(ids.length)label.setAttribute('aria-hidden','true');
  const wrapper=make('span',{className:'number-with-label'});
  wrapper.append(span,label);
  return wrapper;
}
function numericFraction(numerator,denominator) {
  const fraction=make('span',{className:'numeric-fraction'});
  for(const contents of [numerator,denominator]){
    const part=make('span');
    part.append(...[contents].flat());
    fraction.append(part);
  }
  return fraction;
}
function showNumeric(kind,data,line,k,step) {
  document.querySelectorAll('.source-highlight').forEach(node=>node.classList.remove('source-highlight'));
  const host=el('numeric');
  host.replaceChildren();
  const help=el('symbol-help');
  const note=el('numeric-note');
  const breakdown=el('factor-breakdown');
  breakdown.hidden=true; breakdown.replaceChildren();
  help.hidden=true; help.textContent='';
  note.hidden=true; note.textContent='';
  const definitions={
    '[H⁺]':'концентрация ионов водорода в мкмоль/л.',
    C:'общая концентрация кислоты. В балансе зарядов подставляется в мкмоль/л, итоговый обратный расчёт показан в моль/л.',
    N:'числитель. Складывает вклады форм кислоты с коэффициентами 1, 2 и 3 по числу отданных протонов.',
    D:'знаменатель. Сумма вкладов всех форм кислоты, включая недиссоциированную.',
    'K₁':'константа диссоциации первой ступени в мкмоль/л.',
    'K₂':'константа диссоциации второй ступени в мкмоль/л.',
    'K₃':'константа диссоциации третьей ступени в мкмоль/л.',
    'Kₐ':'константа диссоциации соответствующей ступени в мкмоль/л.',
    'pKₐ':'отрицательный десятичный логарифм константы диссоциации.',
    'κ':'электропроводность, измеренная прибором, в мкСм/см.',
    'κ₀':'электропроводность воды или холостой пробы в мкСм/см.',
    'κиспр':'электропроводность после вычитания измеренного фона воды κ₀.',
    'κтеор':'расчётная электропроводность от всех ионов. С автопротолизом выражена как прирост относительно чистой воды.',
    'Kᵥ':'ионное произведение воды при 25 °C. В шкале мкмоль/л равно 0,01 (мкмоль/л)².',
    '[OH⁻]':'концентрация гидроксид-ионов. При включённом автопротолизе равна Kᵥ/[H⁺].',
    'Σλᵢcᵢ':'сумма вкладов ионов. При λ в См·см²/моль и c в мкмоль/л делим на 1000 для результата в мкСм/см.',
    a:'наклон калибровочной прямой.',
    b:'свободный член прямой. Это значение по оси Y при [H⁺] = 0.',
    'R²':'доля разброса показаний, объяснённая линейной регрессией.',
    SSE:'сумма квадратов отклонений показаний от тренда.',
    SST:'сумма квадратов отклонений показаний от среднего.',
  };
  function symbol(label) {
    const button=make('button',{type:'button',className:'symbol-button',textContent:label});
    button.setAttribute('aria-pressed','false');
    button.setAttribute('aria-label','Что означает '+label);
    button.setAttribute('aria-controls','symbol-help');
    button.onclick=()=>{
      host.querySelectorAll('.symbol-button[aria-pressed="true"]').forEach(node=>node.setAttribute('aria-pressed','false'));
      button.setAttribute('aria-pressed','true');
      help.textContent=label+' — '+definitions[label];
      help.hidden=false;
    };
    return button;
  }
  function setNote(text) {note.textContent=text;note.hidden=false;}
  function row(...parts) {
    const p=make('p',{className:'numeric-line'});
    p.append(...parts);
    host.append(p);
  }
  const i=data.i;
  const kSources=records[active].pka.flatMap((_,j)=>['ka-'+j,'pka-'+j]);
  const xSources=model.points.map(p=>(model.lambdas?'theory-':'h-')+p.i);
  const ySources=model.points.map(p=>'reading-'+p.i);
  const fitSources=[...xSources,...ySources,...(model.hasBlank?['blank']:[])];
  const token=numberToken;
  const frac=numericFraction;
  const fixed=(value,label)=>token(value,'fixed',[],label);
  const scale=()=>fixed('10⁶','перевод единиц');
  const divideConductivity=()=>fixed('1000','перевод в мкСм/см');
  const kw=()=>token(WATER.kw,'k',['water-mode'],'Kᵥ при 25 °C');
  const factorLabels=['K₁','K₂','K₃'];
  const stepSources=stage=>stage<records[active].pka.length?['ka-'+stage,'pka-'+stage]:[];
  function openFactors(levels) {
    const label=factorLabels.slice(0,levels).join(' × ');
    const heading=make('h4',{textContent:levels===1?'Как получено K₁':'Как перемножены константы'});
    const expression=make('p',{className:'numeric-line'});
    if(levels===1){
      const pka=parseNumber(records[active].pka[0]);
      expression.append('K₁ = 10^(',fixed(6,'перевод в мкмоль/л'),' − ',
        token(pka,'k',['pka-0'],'pKₐ₁'),') ≈ ',token(k[0],'result',['ka-0'],'K₁, мкмоль/л'));
    }else{
      expression.append(label+' ≈ ');
      for(let j=0;j<levels;j++){
        if(j)expression.append(' × ');
        expression.append(token(k[j],'k',stepSources(j),factorLabels[j]+(stepSources(j).length?'':' · ступени нет')));
      }
      expression.append(' ≈ ',token(k.slice(0,levels).reduce((product,v)=>product*v,1),'result',[], 'Произведение'));
    }
    breakdown.replaceChildren(heading,expression);
    breakdown.hidden=false;
  }
  function kProduct(levels){
    const value=k.slice(0,levels).reduce((product,v)=>product*v,1);
    return token(value,'k',kSources.slice(0,2*levels),factorLabels.slice(0,levels).join(' × '),false,()=>openFactors(levels));
  }
  const concentration=value=>token(value,'input',['conc-'+i],'C стандарта '+(i+1),value<0.001);
  const hValue=()=>token(data.h,'result',['h-'+i],'[H⁺] стандарта '+(i+1));
  const reading=(raw,probe=false)=>token(raw,'signal',[probe?'probe':'reading-'+i],probe?'Показание пробы':'Показание стандарта '+(i+1));
  const water=()=>token(model.blank,'signal',['blank'],'Фон воды κ₀');
  const a=()=>token(line.a,'fit',['slope-value'],'Наклон a');
  const b=()=>token(line.b,'fit',['intercept-value'],'Свободный член b');
  if(kind==='theory'&&data.theory) {
    const t=data.theory;
    const lambdaH=token(WATER.h,'k',['water-mode'],'λ° H⁺');
    const parts=species(data.h,k);
    row('κ H⁺ = ',lambdaH,' × ',hValue(),' / ',divideConductivity(),' ≈ ',
        token(t.contributions[0],'result',['theory-'+i],'Вклад H⁺'),' мкСм/см');
    model.entries.forEach(([name],j)=>{
      row('['+name+'] = ',token(data.c*1e6,'input',['conc-'+i],'C стандарта, мкмоль/л'),' × ',
          token(parts.fractions[j+1],'k',[...kSources,'h-'+i],'Доля '+name),
          ' ≈ ',token(t.concentrations[j],'fit',['h-'+i,'conc-'+i],'Концентрация '+name),' мкмоль/л');
      row('κ '+name+' = ',token(model.lambdas[j],'k',['lambda-'+j],'λ° '+name),' × ',
          token(t.concentrations[j],'fit',['h-'+i,'conc-'+i],'Концентрация '+name),' / ',divideConductivity(),' ≈ ',
          token(t.contributions[j+1],'result',['theory-'+i],'Вклад '+name),' мкСм/см');
    });
    if(model.water)row('κ OH⁻ = ',token(WATER.oh,'k',['water-mode'],'λ° OH⁻'),' × ',
      token(t.oh,'fit',['water-mode','h-'+i],'[OH⁻]'),' / ',divideConductivity(),' ≈ ',
      token(t.contributions.at(-1),'result',['theory-'+i],'Вклад OH⁻'),' мкСм/см');
    const terms=t.contributions.flatMap((value,j)=>[
      ...(j?[' + ']:[]),
      token(value,'fit',['theory-'+i],j===0?'H⁺':model.entries[j-1]?.[0] || 'OH⁻'),
    ]);
    const baseline=model.water?[' − ',token(t.pure,'signal',['water-mode'],'Расчётная чистая вода')]:[];
    row(symbol('κтеор'),' = ',...terms,...baseline,' ≈ ',token(t.value,'result',['theory-'+i],'Расчётная κ'),' мкСм/см');
  }else if(kind==='num'||kind==='den') {
    const context=chargeContext(data,line,k,drillOrigin?.kind);
    if(context){
      const h=()=>token(context.h,'fit',context.ids,'[H⁺] для расчёта');
      const k1=()=>kProduct(1), k12=()=>kProduct(2), k123=()=>kProduct(3);
      if(kind==='num'){
        row(symbol('N'),' = ',k1(),' × ',h(),'² + ',fixed(2,'два протона'),' × ',k12(),' × ',h(),' + ',fixed(3,'три протона'),' × ',k123(),
            ' ≈ ',token(context.parts.num,'result',[],'Числитель N',true));
      }else{
        row(symbol('D'),' = ',h(),'³ + ',k1(),' × ',h(),'² + ',k12(),' × ',h(),' + ',k123(),
            ' ≈ ',token(context.parts.den,'result',[],'Знаменатель D',true));
      }
    }
  }else if(kind==='ka' && k) {
    row(symbol('Kₐ'),' = 10^(',fixed(6,'перевод в мкмоль/л'),' − ',token(parseNumber(records[active].pka[step]),'k',['pka-'+step],'Значение pKₐ'),') = ',
        token(k[step],'result',['ka-'+step],'Пересчитанное значение Kₐ'),' мкмоль/л');
  }else if(kind==='conc' && data.c!==null) {
    row(symbol('C'),' = ',concentration(data.c),' × ',scale(),' = ',token(data.c*1e6,'result',['conc-'+i],'C в мкмоль/л'),' мкмоль/л');
  }else if((kind==='reading'||kind==='blank') && data.raw!==null) {
    if(model.hasBlank) row(symbol('κиспр'),' = ',symbol('κ'),' − ',symbol('κ₀'),' = ',reading(data.measured),' − ',water(),' = ',token(data.y,'result',['reading-'+i,'blank'],'Показание стандарта за вычетом фона'),' мкСм/см');
    else row(symbol('κ'),' = ',reading(data.measured),' мкСм/см. Фон не задан.');
  }else if((kind==='h'||kind==='balance') && data.h!==null && k) {
    const c=data.c*1e6;
    if(model.water||kind==='balance'){
      const values=species(data.h,k);
      const deps=[...kSources,'h-'+i];
      const oh=model.water?WATER.kw/data.h:0;
      if(model.water)row(symbol('[OH⁻]'),' = ',frac(symbol('Kᵥ'),symbol('[H⁺]')),' ≈ ',frac(kw(),hValue()),
        ' ≈ ',token(oh,'result',['water-mode','h-'+i],'[OH⁻] воды'),' мкмоль/л');
      row(symbol('[H⁺]'),' = ',...(model.water?[symbol('[OH⁻]'),' + ']:[]),
          symbol('C'),' × ',frac(symbol('N'),symbol('D')),' ≈ ',
          ...(model.water?[token(oh,'fit',['water-mode','h-'+i],'[OH⁻] воды'),' + ']:[]),
          concentration(c),' × ',frac(token(values.num,'k',deps,'Числитель N',true,'num'),
            token(values.den,'k',deps,'Знаменатель D',true,'den')),' ≈ ',hValue(),' мкмоль/л');
    }else if(k[1]===0){
      const k1=()=>kProduct(1);
      row(symbol('[H⁺]'),' = ',frac(['√(',symbol('K₁'),'² + 4 × ',symbol('K₁'),' × ',symbol('C'),') − ',symbol('K₁')],'2'),
          ' ≈ ',frac(['√(',k1(),'² + ',fixed(4,'коэффициент корня'),' × ',k1(),' × ',concentration(c),') − ',k1()],fixed(2,'делитель')),
          ' ≈ ',hValue(),' мкмоль/л');
    }else{
      const values=species(data.h,k);
      const deps=[...kSources,'h-'+i];
      row(symbol('[H⁺]'),' = ',symbol('C'),' × ',frac(symbol('N'),symbol('D')),' ≈ ',concentration(c),' × ',
          frac(token(values.num,'k',deps,'Числитель N',true,'num'),
            token(values.den,'k',deps,'Знаменатель D',true,'den')),' ≈ ',hValue(),' мкмоль/л');
    }
  }else if(kind==='fit' && line) {
    row(symbol('a'),' = ',token(line.sxy,'input',fitSources,'Сумма по точкам калибровки Sxy',true),' / ',
        token(line.sxx,'input',xSources,'Сумма по оси X · Sxx',true),' ≈ ',a(),' · ',symbol('b'),' = ',
        token(line.yMean,'signal',ySources,'Среднее показание стандартов'),' − ',a(),' × ',
        token(line.xMean,'input',xSources,'Среднее по оси X'),' ≈ ',b());
    row(symbol(model.hasBlank?'κиспр':'κ'),' ≈ ',a(),' × ',symbol(model.lambdas?'κтеор':'[H⁺]'),' + (',b(),')');
  }else if(kind==='r2' && line && line.r2!==null) {
    row(symbol('R²'),' = 1 − ',symbol('SSE'),' / ',symbol('SST'),' ≈ ',fixed(1,'база R²'),' − ',token(line.sse,'input',fitSources,'Сумма квадратов отклонений от тренда',true),' / ',
        token(line.sst,'input',ySources,'Сумма квадратов отклонений показаний от среднего',true),' ≈ ',
        token(line.r2,'result',['metric'],'R² на графике'));
  }else if(kind==='back'||kind==='probe') {
    const probe=kind==='probe';
    const measured=probe?instrumentReading(parseNumber(probeText)):data.measured;
    const result=measured!==null&&measured>=0?invert(measured-model.blank,line,k,model.water,model.lambdas):null;
    if(result){
      const source=probe?'probe':'reading-'+i;
      const hSources=[source,'slope-value','intercept-value',...(model.hasBlank?['blank']:[])];
      const h=()=>token(result.h,'fit',hSources,'[H⁺] пробы');
      const adjusted=()=>token(measured-model.blank,'signal',[source,...(model.hasBlank?['blank']:[])],'Показание с поправкой на фон');
      if(model.hasBlank) row(symbol('κиспр'),' = ',symbol('κ'),' − ',symbol('κ₀'),' = ',reading(measured,probe),' − ',water(),' = ',adjusted(),' мкСм/см');
      const xSymbol=model.lambdas?'κтеор':'[H⁺]';
      const xToken=()=>model.lambdas?token(result.x,'fit',hSources,'κтеор по тренду'):h();
      row(symbol(xSymbol),' = (',symbol(model.hasBlank?'κиспр':'κ'),' − ',symbol('b'),') / ',symbol('a'),' ≈ (',adjusted(),' − (',b(),')) / ',a(),' ≈ ',xToken(),model.lambdas?' мкСм/см':' мкмоль/л');
      if(model.lambdas)row(symbol('Σλᵢcᵢ'),model.water?' минус вклад чистой воды':'',' = ',xToken(),
        ' мкСм/см → ',symbol('[H⁺]'),' ≈ ',h(),' мкмоль/л');
      const oh=()=>token(result.oh,'fit',['water-mode',...hSources],'[OH⁻] воды');
      row(symbol('C'),' = ',...(model.water?['(',symbol('[H⁺]'),' − ',symbol('[OH⁻]'),')']:[symbol('[H⁺]')]),
          ' × ',frac(symbol('D'),symbol('N')),' ÷ 10⁶ ≈ ',...(model.water?['(',h(),' − ',oh(),')']:[h()]),' × ',frac(token(result.den,'k',[...hSources,...kSources],'Знаменатель D',true,'den'),
            token(result.num,'k',[...hSources,...kSources],'Числитель N',true,'num')),' ÷ ',scale(),' ≈ ',
          token(result.c,'result',probe?[]:['back-'+i],'Концентрация по тренду',true),' моль/л');
      if(!probe&&data.c>0)setNote('Отклонение от стандарта '+fmt(100*(result.c/data.c-1),4)+'%.');
      if(probe && model.points.length && (measured<Math.min(...model.points.map(p=>p.measured))||measured>Math.max(...model.points.map(p=>p.measured)))) setNote('Показание вне диапазона стандартов.');
    }else if(probe&&!probeText)row('Введи показание пробы.');
  }
  if(!host.childElementCount) row('Числа появятся, когда будут корректные исходные данные.');
}
el('formula-back').onclick=()=>{
  if(drillOrigin){
    const previous=drillOrigin;
    drillOrigin=null;
    select(previous.kind,previous.row,previous.step);
    el('mode').focus({preventScroll:true});
  }
};
el('mode').onchange=event=>select(event.target.value);
el('water-mode').onchange=refresh;
el('ion-mode').onchange=refresh;
el('model-toggle').onchange=event=>{
  advanced=event.currentTarget.checked;
  el('physics').hidden=!advanced;
  el('physics-explanation').hidden=!advanced;
  el('water-mode').checked=advanced;
  el('ion-mode').checked=advanced;
  el('ion-editor').open=false;
  if(drillOrigin){
    selected=drillOrigin;
    drillOrigin=null;
    paintSelection();
  }
  refresh();
};

el('instrument-toggle').checked=instrumentCorrectionCancelled;
el('instrument-toggle').onchange=event=>{
  instrumentCorrectionCancelled=event.currentTarget.checked;
  try { localStorage.setItem(instrumentKey,String(instrumentCorrectionCancelled)); } catch {}
  el('probe').value=displayReading(probeText);
  mount();
};
el('probe').oninput=event=>{probeText=inputReading(event.target.value); showFormula();};
el('probe').onblur=event=>{event.target.value=displayReading(probeText);};
el('metric').onclick=()=>select('r2');
el('trend').onclick=()=>select('fit');
window.addEventListener('scroll',hidePointTooltip,true);
window.addEventListener('resize',hidePointTooltip);
document.addEventListener('pointerdown',event=>{
  if(!event.target.closest('.point,.chart-dot'))hidePointTooltip();
});
document.addEventListener('keydown',event=>{
  if(event.key==='Escape')hidePointTooltip();
});
el('save-sheet').onclick=()=>{
  el('save-status').textContent=persist()?'Сохранено в этом браузере.':'Браузер не разрешил сохранить данные.';
};
function exportChoice(value,text,checked,quality) {
  const label=make('label');
  if(quality)label.dataset.quality=quality;
  label.append(make('input',{type:'checkbox',value,checked}),text);
  return label;
}
const exportColors=Object.entries(diagnosticLabels).map(([quality,label])=>
  exportChoice(quality,label,['warning','influence','repeat'].includes(quality),quality));
el('export-colors').append(...exportColors);
const exportAcids=el('export-acids');
const checkedValues=container=>[...container.querySelectorAll('input:checked')].map(input=>input.value);
el('export-errors').onclick=()=>{
  const previous=new Map([...exportAcids.querySelectorAll('input')].map(input=>[input.value,input.checked]));
  exportAcids.replaceChildren(...names.slice().sort((a,b)=>a.localeCompare(b,'ru'))
    .map(name=>exportChoice(name,acidTitle(name),previous.get(name)??true)));
  el('export-notice').textContent='';
  el('export-dialog').showModal();
};
el('export-cancel').onclick=()=>el('export-dialog').close();
for(const [id,checked] of [['export-all',true],['export-none',false]])
  el(id).onclick=()=>exportAcids.querySelectorAll('input').forEach(input=>{input.checked=checked;});
function calibrationExportRows(acids,colors,includePlanned) {
  return acids.flatMap(name=>{
    const {rows,points,line,blank}=compute(name);
    const selected=points.filter(point=>colors.has(point.diagnostics.quality));
    if(includePlanned&&line){
      const concentrations=new Set(records[name].hexBnPlanned);
      selected.push(...rows.filter(row=>row.raw===null&&row.x!==null&&
        String(records[name].rows[row.i][1]).trim()===''&&concentrations.has(row.c)));
    }
    return selected.sort((a,b)=>b.c-a.c).map(point=>{
      const predicted=line?line.a*point.x+line.b+blank:null;
      return [acidTitle(name),point.c,Number.isFinite(predicted)?roundReading(predicted):null,point.measured];
    });
  });
}
el('export-form').onsubmit=event=>{
  event.preventDefault();
  const acids=checkedValues(exportAcids), colors=new Set(checkedValues(el('export-colors')));
  const notice=el('export-notice');
  const includePlanned=el('export-planned').checked;
  if(!acids.length||(!colors.size&&!includePlanned)){
    notice.textContent='Выбери хотя бы одну кислоту и отметь цвета точек или плановые стандарты.';
    return;
  }
  const rows=calibrationExportRows(acids,colors,includePlanned);
  if(!rows.length){
    notice.textContent='Для выбранных кислот подходящих строк нет.';
    return;
  }
  const file=calibrationWorkbook(rows);
  const url=URL.createObjectURL(file);
  const link=make('a',{href:url,download:'calibration-errors.xlsx'});
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(()=>URL.revokeObjectURL(url),1000);
  notice.textContent='Выгружено точек: '+rows.length+'.';
};
el('add-row').onclick=()=>{
  const rows=records[active].rows;
  rows.push(['','']);
  persist();
  showInputs();
  refresh();
  el('conc-'+(rows.length-1)).focus();
};
el('add-acid').onclick=()=>{
  el('acid-form').reset();
  el('acid-name').setCustomValidity('');
  el('acid-dialog').showModal();
};
el('cancel-acid').onclick=()=>el('acid-dialog').close();
el('acid-name').oninput=()=>el('acid-name').setCustomValidity('');
el('acid-form').onsubmit=event=>{
  event.preventDefault();
  const input=el('acid-name');
  const name=input.value.trim();
  const duplicate=names.some(existing=>
    existing.toLocaleLowerCase('ru')===name.toLocaleLowerCase('ru') ||
    acidTitle(existing).toLocaleLowerCase('ru')===name.toLocaleLowerCase('ru'));
  if(!name||duplicate){
    input.setCustomValidity(name?'Кислота с таким названием уже есть.':'Введи название кислоты.');
    input.reportValidity();
    return;
  }
  records[name]=newAcid(Number(el('acid-steps').value));
  names.push(name);
  active=name;
  selected={kind:'h',row:0,step:0};
  drillOrigin=null;
  persist();
  el('acid-dialog').close();
  mount();
  focusAcid(name);
  el('pka-0').focus();
};
function saveTabColor(color) {
  if(color)records[colorAcid].color=color;
  else delete records[colorAcid].color;
  persist();
  el('color-dialog').close();
  refresh();
  focusAcid(colorAcid,'.tab-menu');
}
el('color-form').onsubmit=event=>{
  event.preventDefault();
  saveTabColor(el('tab-color').value);
};
el('cancel-color').onclick=()=>el('color-dialog').close();
el('reset-color').onclick=()=>saveTabColor('');
el('delete-acid').onclick=()=>{
  const name=colorAcid;
  if(Object.hasOwn(seed,name)||!confirm('Удалить кислоту «'+name+'» и все её данные?'))return;
  const index=names.indexOf(name);
  names.splice(index,1);
  delete records[name];
  persist();
  if(active===name){
    active=names[Math.min(index,names.length-1)];
    selected={kind:'h',row:0,step:0};
    drillOrigin=null;
  }
  el('color-dialog').close();
  mount();
  focusAcid(names[Math.min(index,names.length-1)]);
};
el('reset').onclick=()=>{
  const original=Object.hasOwn(seed,active);
  const message=original?'Восстановить исходные точки и pKa на листе':'Очистить введённые данные на листе';
  if(confirm(message+' «'+active+'»?')){
    const color=tabColor(active);
    records[active]=original?clone(seed[active]):newAcid(records[active].pka.length);
    if(color)records[active].color=color;
    persist();
    selected={kind:'h',row:0,step:0};
    drillOrigin=null;
    mount();
  }
};
const fractionParts = [
  {key:'h',description:'h = [H⁺] — равновесная концентрация ионов водорода. Для обратного расчёта подставляем h, полученную по калибровке, в мкмоль/л.'},
  {key:'c',description:'C — аналитическая концентрация кислоты, которую определяем по [H⁺]. Для стандарта её можно сравнить со значением в столбце A.'},

  {key:'n',description:'N — знаменатель в формуле для C. Каждый вклад умножен на число протонов, отданных соответствующей формой кислоты: 1, 2 или 3.'},
  {key:'d',description:'D — числитель в формуле для C. Это сумма вкладов всех форм кислоты, включая недиссоциированную. Отношение N/D показывает среднее число отданных протонов на молекулу кислоты.'},
];
let chosenFraction = null;
function showFraction(active) {
  fractionParts.forEach(entry=>{
    const highlighted=String(entry===active);
    entry.parts.forEach(part=>{part.dataset.highlight=highlighted;});
    entry.button.dataset.highlight=highlighted;
    entry.button.setAttribute('aria-pressed',String(entry===chosenFraction));
  });
  el('fraction-explanation').textContent=active?.description || 'Выбери обозначение, чтобы увидеть его смысл.';
}
for(const entry of fractionParts){
  entry.parts=[...document.querySelectorAll('[data-symbol='+entry.key+']')];
  entry.button=el('explain-'+entry.key);
  for(const target of [...entry.parts,entry.button]){
    target.addEventListener('pointerenter',()=>showFraction(entry));
    target.addEventListener('pointerleave',()=>showFraction(chosenFraction));
    target.addEventListener('click',()=>{
      chosenFraction=entry;
      showFraction(entry);
    });
  }
}
appendHexBnStandards();
if(location.hash==='#hexbn-calibrations'&&window.opener){
  const hexOrigin=location.origin;
  let savedMode={};
  try { savedMode=JSON.parse(localStorage.getItem('ibu-bn-model-v1')||'{}')||{}; } catch {}
  window.opener.postMessage({
    type:'hexbn-calibrations-v1',records,water:savedMode.water===true,ions:savedMode.ions===true
  },hexOrigin);
  window.close();
} else {
  mount();
}
