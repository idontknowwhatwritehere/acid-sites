const origin='https://ibu-bn-acid-calibration.nocanbewant.chatgpt.site';
const key='hexbn-experiment-v1';
const base={vd:80,vr:30,cd:0.01,cr:0,area:10.57309,extrapolate:true};
const fields=[['vd','Объём донора, мл'],['vr','Объём приёмника, мл'],['cd','C донора при t = 0, М'],['cr','C приёмника при t = 0, М'],['area','Площадь мембраны, см²']];
const el=id=>document.getElementById(id);
const fmt=(x,d=5)=>Number.isFinite(x)?new Intl.NumberFormat('ru',{maximumSignificantDigits:d}).format(x):'—';
const sci=x=>Number.isFinite(x)?x.toExponential(4).replace('.',','):'—';
const val=x=>String(x??'').replace('.',',');
function node(tag,label){const e=document.createElement(tag);if(label!==undefined)e.textContent=label;return e;}
let state={edits:{},params:{},sync:null};
try {const saved=JSON.parse(localStorage.getItem(key));if(saved)state=saved;} catch {}
state.edits=Object.assign(Object.create(null),state.edits);
state.params=Object.assign(Object.create(null),state.params);
let acids=[],acid=state.acid,series=state.series??0,selected=0,
  view=state.view==='calibration'?'calibration':'permeability',syncWindow,syncTimer;
function save(){
  Object.assign(state,{acid,series,view});
  try {localStorage.setItem(key,JSON.stringify(state));}
  catch {el('sync-status').textContent='Браузер не смог сохранить изменения.';}
}
const records=()=>state.sync?.records??defaults;
const runs=()=>Object.hasOwn(experiments,acid)?experiments[acid]:[{rows:[],start:null,end:null}];
const params=()=>({...base,...state.params[acid]});
const rawRows=()=>state.edits[acid]?.[series]??runs()[series].rows;
function editRows(){const entries=state.edits[acid]??(state.edits[acid]=[]);return entries[series]??(entries[series]=structuredClone(rawRows()));}
function selectOptions(select,options,value){
  select.replaceChildren(...options.map(([key,text])=>new Option(text,key)));
  select.value=value;
}
function mount(){
  acids=[...new Set([...Object.keys(experiments),...Object.keys(records())])].sort((a,b)=>a.localeCompare(b,'ru'));
  if(!acids.includes(acid))acid=acids[0];
  if(!runs()[series])series=0;
  selectOptions(el('acid'),acids.map(n=>[n,n]),acid);
  selectOptions(el('series'),runs().map((r,i)=>[i,'Серия '+(i+1)]),series);
  const cal=view==='calibration';
  el('series-control').hidden=cal;
  el('parameters').hidden=cal;
  el('actions').hidden=cal;
  el('calibration-view').setAttribute('aria-pressed',String(cal));
  el('permeability-view').setAttribute('aria-pressed',String(!cal));
  const ref=runs()[series];
  const edited=Array.isArray(state.edits[acid]?.[series]);
  el('provenance').textContent=cal?'i‑BuBn':edited?'Изменено в браузере':ref.start?'Кислоты1.xlsx · A'+ref.start+':B'+ref.end:'Нет исходных измерений';
  el('sheet-title').textContent=cal?'Калибровочные стандарты':'Измерения HexBn';
  el('extrapolate').checked=params().extrapolate;
  el('parameters').replaceChildren(...fields.map(([id,label])=>{
    const field=node('label',label), input=node('input');
    input.value=val(params()[id]);input.inputMode='decimal';input.id='param-'+id;
    input.oninput=()=>{(state.params[acid]??={})[id]=input.value;save();refresh();};
    field.append(input);return field;
  }));
  buildTable();
  refresh();
}
function calibration(){
  const record=Object.hasOwn(records(),acid)?records()[acid]:null;
  const k=record?constants(record.pka):null;
  const blank=record?Math.max(0,parseNumber(record.blank??1.7)??0):0;
  const water=state.sync?.water===true;
  const entries=Object.hasOwn(ionData.ions,acid)?ionData.ions[acid]:
    (record?.pka??[]).map(()=>[null,null]);
  const values=entries.map((entry,i)=>parseNumber(record?.lambda?.[i]??entry[1]));
  const lambdas=state.sync?.ions&&values.length&&values.every(v=>v>0)?values:null;
  const rows=(record?.rows??[]).map(([concentration,reading],i)=>{
    const c=parseNumber(concentration), raw=parseNumber(reading);
    const h=k&&c>0?hydrogen(c,k,water):null;
    const x=h>0&&lambdas?ionicConductivity(h,c,species(h,k),water,lambdas).value:h;
    return {i,c,raw,h,x,y:raw===null?null:raw-blank};
  });
  const points=rows.filter(p=>Number.isFinite(p.x)&&p.x>0&&p.raw!==null&&p.raw>=0);
  return {rows,points,line:fit(points),blank,k,water,lambdas,
    min:Math.min(...points.map(p=>p.raw)),max:Math.max(...points.map(p=>p.raw))};
}
function permeability(cal){
  const p=params();
  const [vd,vr,cd,cr,area]=fields.map(([id])=>parseNumber(p[id]));
  const valid=vd>0&&vr>0&&area>0&&cr!==null&&cr>=0&&cd>cr;
  const factor=valid?area*(1/vd+1/vr):null;
  const rows=rawRows().map(([time,reading],i)=>{
    const t=parseNumber(time), raw=parseNumber(reading);
    const outside=raw!==null&&cal.points.length>0&&(raw<cal.min||raw>cal.max);
    const back=raw!==null&&raw>=0?invert(raw-cal.blank,cal.line,cal.k,cal.water,cal.lambdas):null;
    const receiver=back?.c??null;
    const donor=valid&&receiver!==null?cd-vr/vd*(receiver-cr):null;
    const delta=donor!==null?donor-receiver:null;
    const fraction=delta!==null?delta/(cd-cr):null;
    const z=fraction>0?-Math.log(fraction):null;
    const eligible=valid&&t!==null&&t>=0&&Number.isFinite(z)&&(!outside||p.extrapolate);
    const transferred=valid&&receiver!==null?vr/1000*(receiver-cr):null;
    return {i,t,raw,outside,back,receiver,donor,delta,z,eligible,transferred};
  });
  const points=rows.filter(r=>r.eligible).map(r=>({...r,x:r.t*60,y:r.z}));
  const n=points.length, sxx=points.reduce((sum,r)=>sum+r.x**2,0);
  const beta=n>=2&&sxx>0?points.reduce((sum,r)=>sum+r.x*r.y,0)/sxx:null;
  const sse=beta===null?null:points.reduce((sum,r)=>sum+(r.y-beta*r.x)**2,0);
  const mean=n?points.reduce((sum,r)=>sum+r.y,0)/n:0;
  const sst=points.reduce((sum,r)=>sum+(r.y-mean)**2,0);
  const r2=sse!==null&&sst>0?1-sse/sst:null;
  const result=beta!==null&&beta>=0?beta/factor:null;
  const se=result!==null?Math.sqrt(sse/(n-1)/sxx)/factor:null;
  return {rows,points,beta,r2,result,se,sse,sxx,rmse:sse===null?null:Math.sqrt(sse/n),valid,factor,vd,vr,cd,cr,area};
}
function buildTable(){
  const cal=view==='calibration', table=el('sheet');
  const header=node('tr');
  for(const label of (cal?['№','C, М','κ, мкСм/см','X модели','κ − κ₀']:['№','t, мин','κ, мкСм/см','C приёмника, М','C донора, М','Статус'])){
    const th=node('th',label);th.scope='col';header.append(th);
  }
  table.tHead.replaceChildren(header);
  const rows=cal?calibration().rows:rawRows();
  table.tBodies[0].replaceChildren(...rows.map((r,i)=>{
    const tr=node('tr'),th=node('th');th.scope='row';
    const pick=node('button',i+1);pick.type='button';pick.title='Показать подстановку';pick.onclick=()=>choose(i);
    th.append(pick);tr.append(th);
    if(!cal){
      for(const j of [0,1]){
        const td=node('td'),input=node('input');
        input.value=val(r[j]);input.inputMode='decimal';
        input.setAttribute('aria-label',(j===0?'Время':'Проводимость')+' · строка '+(i+1));
        input.onfocus=()=>choose(i);
        input.oninput=()=>{editRows()[i][j]=input.value;save();refresh();};
        td.append(input);tr.append(td);
      }
    }
    for(let j=0;j<(cal?4:3);j++){const td=node('td');td.className='output';tr.append(td);}
    if(!cal){
      const remove=node('button','×');remove.type='button';remove.className='remove-row';
      remove.setAttribute('aria-label','Удалить строку '+(i+1));remove.onclick=()=>{
        if(rawRows()[i].some(v=>String(v).trim())&&!confirm('Удалить строку '+(i+1)+'?'))return;
        editRows().splice(i,1);selected=Math.max(0,Math.min(i,rawRows().length-1));save();mount();
      };
      th.append(remove);
    }
    return tr;
  }));
}
let current;
function choose(i){selected=i;paintDetail();}
function paintDetail(){
  if(!current)return;
  const {cal,run}=current, isCal=view==='calibration';
  const r=(isCal?cal.rows:run.rows)[selected];
  for(const [i,tr] of [...el('sheet').tBodies[0].rows].entries())tr.dataset.selected=String(i===selected);
  for(const point of document.querySelectorAll('#plot [data-row], #other-plots [data-row]'))point.setAttribute('aria-pressed',String(+point.dataset.row===selected));
  renderFormulas(cal,run,r,isCal);
  el('detail-title').textContent=r?'Подстановка · строка '+(selected+1):'Подстановка для точки';
  if(!r){el('detail').textContent='Выбери точку графика или номер строки.';return;}
  const a=fmt(cal.line?.a),b=fmt(cal.line?.b);
  const calibrationText='κ − κ₀ = '+fmt(r.raw)+' − '+fmt(cal.blank)+' = '+fmt(r.raw===null?null:r.raw-cal.blank)+' мкСм/см';
  const parts=isCal?[
    'C стандарта = '+sci(r.c)+' М',calibrationText,
    (cal.lambdas?'κ расчётная':'[H⁺]')+' = '+fmt(r.x)+(cal.lambdas?' мкСм/см':' мкмоль/л')
  ]:[
    't = '+fmt(r.t)+' мин; '+calibrationText,
    'X = ((κ − κ₀) − b) / a; a = '+a+', b = '+b+'; X = '+fmt(r.back?.x),
    'C приёмника = '+sci(r.receiver)+' М',
    'C донора = '+fmt(run.cd)+' − '+fmt(run.vr)+'/'+fmt(run.vd)+' × ('+sci(r.receiver)+' − '+fmt(run.cr)+') = '+sci(r.donor)+' М',
    'z = −ln[(C донора − C приёмника)/(C донора,0 − C приёмника,0)] = '+fmt(r.z),
    'P = β / [A × (1/Vд + 1/Vп)] = '+sci(run.beta)+' / '+sci(run.factor)+' = '+sci(run.result)+' см/с'
  ];
  el('detail').textContent=parts.join('\n');
}
function refresh(){
  const isCal=view==='calibration',cal=calibration(),run=permeability(cal);
  current={cal,run};
  el('result').textContent=isCal?'R² '+fmt(cal.line?.r2):'P '+sci(run.result)+' см/с';
  const seriesStats=isCal?[
    ['Точек',cal.points.length],['Наклон a',fmt(cal.line?.a)],['Свободный член b',fmt(cal.line?.b)],['Фон κ₀',fmt(cal.blank)]
  ]:[
    ['В регрессии',run.points.length+' / '+run.rows.length],['Из них экстраполяций',run.points.filter(r=>r.outside).length],
    ['R²',fmt(run.r2)],['RMSE по z',fmt(run.rmse)],['β, с⁻¹',sci(run.beta)],['Ст. ошибка P, см/с',sci(run.se)]
  ];
  el('statistics').replaceChildren(...seriesStats.map(([label,value])=>{const wrap=node('div');wrap.append(node('dt',label),node('dd',value));return wrap;}));
  const messages=[];
  if(!cal.line||!(cal.line.a>0))messages.push('Для обратного расчёта нужна калибровка с положительным наклоном.');
  if(!isCal){
    if(!run.valid)messages.push('Проверь объёмы, площадь и начальные концентрации.');
    if(run.points.length<4)messages.push('Мало пригодных точек. Оценка P чувствительна к отдельным измерениям.');
    if(run.beta<0)messages.push('Отрицательный наклон. Проверь направление переноса и входные данные.');
    if(run.rows.some(r=>r.outside))messages.push('Оранжевые строки выходят за диапазон калибровки '+fmt(cal.min)+'–'+fmt(cal.max)+' мкСм/см.');
  }
  el('fit-note').textContent=messages.join(' ');
  (isCal?cal.rows:run.rows).forEach((r,i)=>{
    const tr=el('sheet').tBodies[0].rows[i];
    const cells=tr.querySelectorAll('.output');
    tr.dataset.outside=String(!isCal&&r.outside);
    const readings=isCal?[sci(r.c),fmt(r.raw),fmt(r.x),fmt(r.y)]:[
      sci(r.receiver),sci(r.donor),r.eligible?(r.outside?'Экстраполяция':'В расчёте'):
      r.outside&&!params().extrapolate?'Вне диапазона':r.raw===null||r.t===null?'Нет пары t, κ':
      r.receiver===null?'C не рассчитана':r.donor!==null&&r.z===null?'Градиент ≤ 0':'Исключена'
    ];
    readings.forEach((value,j)=>cells[j].textContent=value);
  });
  showHexPlot(cal,run,isCal);
  el('sync-status').textContent=state.sync?'Калибровки получены '+new Date(state.sync.at).toLocaleString('ru-RU')+'.':'Исходные калибровки i‑BuBn. Чтобы учесть правки из своего браузера, подтяни калибровки.';
  el('method').innerHTML=isCal?
    '<h3>Калибровочная прямая</h3><p class="equation">κ − κ₀ = aX + b</p><p>X — '+(cal.lambdas?'расчётная проводимость по сумме вкладов ионов.':'концентрация H⁺ из баланса зарядов с учётом ступеней диссоциации.')+' a — наклон, b — свободный член, κ₀ — фон воды. Для обратного расчёта сначала находят X = (κ − κ₀ − b) / a, затем концентрацию кислоты C из того же равновесия.</p>':
    '<h3>Модель двух камер с конечными объёмами</h3><p class="equation">Cд = Cд,0 − (Vп/Vд)(Cп − Cп,0)</p><p>Калибровка переводит измеренную проводимость κ в концентрацию приёмника Cп. По балансу вещества вычисляем Cд в отдающей фазе. Площадь по умолчанию 10,57309 см², если A1 книги задана в м².</p><p class="equation">z = −ln[(Cд − Cп)/(Cд,0 − Cп,0)] = βt</p><p>t в секундах. β — наклон прямой через начало координат, β = Σ(tz)/Σ(t²). Такой расчёт привязывает начальный градиент к введённым начальным концентрациям.</p><p class="equation">P = β / [A(1/Vд + 1/Vп)]</p><p>Vд и Vп в мл (см³), A в см², P в см/с. Модель предполагает хорошо перемешанные камеры, постоянные объёмы и отсутствие потерь вещества. Стандартная ошибка P отражает разброс относительно тренда. Неопределённость калибровки в неё не входит. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8951150/" target="_blank" rel="noopener">Источник модели</a>.</p>';
  paintDetail();
}
el('acid').onchange=e=>{acid=e.target.value;series=0;selected=0;save();mount();};
el('series').onchange=e=>{series=Number(e.target.value);selected=0;save();mount();};
el('permeability-view').onclick=()=>{view='permeability';selected=0;save();mount();};
el('calibration-view').onclick=()=>{view='calibration';selected=0;save();mount();};
el('extrapolate').onchange=e=>{(state.params[acid]??={}).extrapolate=e.target.checked;save();refresh();};
el('add').onclick=()=>{const rows=editRows();rows.push(['','']);selected=rows.length-1;save();mount();el('sheet').tBodies[0].rows[selected].querySelector('input').focus();};
el('paste').onclick=()=>{el('paste-form').reset();el('paste-error').textContent='';el('paste-dialog').showModal();};
el('paste-cancel').onclick=()=>el('paste-dialog').close();
el('paste-form').onsubmit=e=>{
  e.preventDefault();
  const rows=el('paste-text').value.trim().split(/\r?\n/).filter(s=>s.trim()).map(s=>s.split(/\t|;/).slice(0,2).map(v=>v.trim()));
  if(rows[0]?.length===2&&rows[0].every(v=>parseNumber(v)===null)&&rows[0].every(v=>/\p{L}/u.test(v)))rows.shift();
  if(!rows.length||rows.some(r=>r.length!==2||r.some(v=>parseNumber(v)===null))){
    el('paste-error').textContent='Нужны два числовых столбца. Проверь время и проводимость во всех строках.';return;
  }
  (state.edits[acid]??=[])[series]=rows;
  selected=0;save();el('paste-dialog').close();mount();
};
el('sync').onclick=()=>{
  clearTimeout(syncTimer);
  syncWindow=window.open(origin+'/#hexbn-calibrations','hex-calibration-sync');
  el('sync-status').textContent=syncWindow?'Ожидаю калибровки из i‑BuBn…':'Разреши всплывающее окно, затем повтори синхронизацию.';
  if(syncWindow)syncTimer=setTimeout(()=>{el('sync-status').textContent='Калибровки не получены. Войди в i‑BuBn в этом браузере и повтори.';},20000);
};
window.addEventListener('message',e=>{
  if(e.origin!==origin||!syncWindow||e.source!==syncWindow||e.data?.type!=='hexbn-calibrations-v1')return;
  const {records:received,water,ions}=e.data;
  if(!received||typeof received!=='object')return;
  clearTimeout(syncTimer);syncWindow=null;
  state.sync={records:received,water:water===true,ions:ions===true,at:Date.now()};
  save();mount();
});
mount();
