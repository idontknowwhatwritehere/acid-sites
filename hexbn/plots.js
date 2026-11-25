const hexPlots=[
  ['conductivity','Электропроводность κ(t)','κ, мкСм/см','Исходные показания прибора до вычитания фона.'],
  ['receiver','Концентрация приёмника Cп(t)','Cп, моль/л','Обратный расчёт по калибровке выбранной кислоты.'],
  ['donor','Концентрация донора Cд(t)','Cд, моль/л','Концентрация в отдающей фазе по балансу вещества.'],
  ['delta','Разность концентраций ΔC(t)','ΔC, моль/л','Разность Cд − Cп.'],
  ['intervalP','Интервальная проницаемость P(t)','P, см/с','Наклон между соседними пригодными замерами. По оси времени — середина интервала.'],
  ['flux','Поток через мембрану J(t)','J, моль/(см²·с)','Средний поток между двумя соседними замерами. Время соответствует середине интервала.'],
  ['linear','Логарифмическая линеаризация z(t)','z','Точки, вошедшие в основной расчёт. Прямая z = βt проходит через начало координат.'],
  ['residual','Остатки регрессии ε(t)','ε = z − βt','Разность между измеренным z и значением тренда.'],
  ['start','Проницаемость P(tнач)','P, см/с','Для каждого фактического начала участка заново оценивается наклон по последующим точкам относительно z и t в начале. Нужны две последующие пригодные точки.'],
  ['fraction','Степень переноса E(t)','E, %','Доля начального количества кислоты в доноре, перешедшая в приёмник.'],
  ['moles','Перенесённое вещество n(t)','n, моль','Прирост количества кислоты в приёмнике относительно начала опыта.'],
  ['fluxDelta','Поток J(ΔC)','J, моль/(см²·с)','Для каждого интервала по X — логарифмическое среднее начальной и конечной ΔC, по Y — средний поток. Линия показывает J = P × ΔC / 1000 для общего P.'],
  ['sse','Сумма квадратов SSE(P)','SSE по z','Сумма (z − P × A × (1/Vд + 1/Vп) × t)². Отмечена оценка P по основной регрессии.'],
  ['cook','Расстояние Кука','Расстояние Кука D','Влияние каждой точки на регрессию через начало координат. Пунктир — ориентир 4/n.']
];
function plotNumber(value){
  return Number.isFinite(value)&&value!==0&&(Math.abs(value)<0.001||Math.abs(value)>=10000)?sci(value):fmt(value,4);
}
function intervalDiagnostics(run){
  const ordered=run.rows.filter(r=>Number.isFinite(r.t)&&r.t>=0).sort((a,b)=>a.t-b.t);
  const intervals=[];
  for(let i=1;i<ordered.length;i++){
    const a=ordered[i-1],b=ordered[i],dt=(b.t-a.t)*60;
    if(!a.eligible||!b.eligible||!(dt>0))continue;
    const change=b.z-a.z;
    const gradient=change===0?a.delta:(a.delta-b.delta)/change;
    const dn=b.transferred-a.transferred;
    intervals.push({i:b.i,x:(a.t+b.t)/2,p:change/(run.factor*dt),dt,dn,dz:change,
      flux:dn/(run.area*dt),gradient,outside:a.outside||b.outside,
      span:fmt(a.t)+'–'+fmt(b.t)+' мин'});
  }
  return intervals;
}
function hexPlotData(run,id){
  const time='Время, мин';
  const samples=key=>run.rows.filter(r=>r.t>=0&&Number.isFinite(r.t)&&Number.isFinite(r[key])&&r.raw>=0)
    .map(r=>({i:r.i,x:r.t,y:r[key],outside:key!=='raw'&&r.outside}));
  const fitPoints=run.points.map(p=>({...p,x:p.t}));
  const interval=intervalDiagnostics(run);
  const basic={conductivity:'raw',receiver:'receiver',donor:'donor',delta:'delta',moles:'transferred'};
  if(Object.hasOwn(basic,id))return {points:samples(basic[id]),xlabel:time};
  if(id==='fraction')return {points:samples('transferred').map(p=>({...p,y:100*p.y/(run.vd/1000*run.cd)})),xlabel:time};
  if(id==='linear')return {points:fitPoints,xlabel:time,line:run.beta===null?null:{a:60*run.beta,b:0}};
  if(id==='residual')return {points:run.beta===null?[]:fitPoints.map(p=>({...p,y:p.z-60*run.beta*p.t})),xlabel:time,line:{a:0,b:0}};
  if(id==='intervalP'||id==='flux'||id==='fluxDelta'){
    return {points:interval.map(p=>({...p,x:id==='fluxDelta'?p.gradient:p.x,y:id==='intervalP'?p.p:p.flux})),
      xlabel:id==='fluxDelta'?'Средняя ΔC, моль/л':time,
      line:id==='fluxDelta'&&run.result!==null?{a:run.result/1000,b:0}:null};
  }
  if(id==='start'){
    const order=fitPoints.slice().sort((a,b)=>a.t-b.t);
    const points=order.flatMap(start=>{
      const later=order.filter(r=>r.t>start.t);
      if(later.length<2)return [];
      const parts=later.map(r=>({t:(r.t-start.t)*60,z:r.z-start.z}));
      const beta=parts.reduce((s,r)=>s+r.t*r.z,0)/parts.reduce((s,r)=>s+r.t*r.t,0);
      return [{...start,y:beta/run.factor,outside:start.outside||later.some(r=>r.outside)}];
    });
    return {points,xlabel:'Начало рабочего участка, мин'};
  }
  if(id==='sse'){
    if(run.result===null)return {points:[],xlabel:'P, см/с'};
    const top=Math.max(run.result*2,run.se*6,1e-12);
    const sse=p=>run.points.reduce((s,r)=>s+(r.z-p*run.factor*r.x)**2,0);
    const curve=Array.from({length:51},(_,i)=>({x:top*i/50,y:sse(top*i/50)}));
    return {points:[{x:run.result,y:sse(run.result),span:'Минимум SSE'}],curve,xlabel:'P, см/с'};
  }
  const n=fitPoints.length,mse=n>1?run.sse/(n-1):0;
  const points=run.beta===null||!(mse>0)?[]:fitPoints.flatMap(p=>{
    const h=(60*p.t)**2/run.sxx, residual=p.z-60*run.beta*p.t;
    return h<1-1e-12?[{...p,x:p.i+1,y:residual**2/mse*h/(1-h)**2}]:[];
  });
  return {points,xlabel:'Номер строки в таблице',line:n?{a:0,b:4/n}:null};
}
function showHexPlot(cal,run,isCal){
  const main=el('plot'),other=el('other-plots');
  el('plot-title').textContent=isCal?'График калибровки':'Проницаемость';
  main.replaceChildren();
  other.replaceChildren();
  other.hidden=isCal;
  if(isCal){
    const graph=node('div');graph.className='plot-wrap';
    main.append(graph);
    renderHexPlot(graph,{points:cal.points,line:cal.line,
      xlabel:cal.lambdas?'κ расчётная, мкСм/см':'[H⁺], мкмоль/л'},'κ − κ₀, мкСм/см');
    return;
  }
  const panel=(id,card=false)=>{
    const [,title,ylabel,note]=hexPlots.find(([value])=>value===id);
    const article=node('article');article.className=card?'card hex-plot-panel':'hex-plot-panel';
    const graph=node('div');graph.className='plot-wrap';
    const description=node('p',note);description.className='hex-foot';
    article.append(node('h3',title),graph,description);
    renderHexPlot(graph,hexPlotData(run,id),ylabel);
    return article;
  };
  main.append(panel('intervalP'),panel('start'));
  const groups=[
    ['Измерения и перенос',['conductivity','receiver','donor','delta','moles','fraction','flux','fluxDelta','linear']],
    ['Статистика регрессии',['residual','sse','cook']]
  ];
  for(const [title,ids] of groups){
    const section=node('section');section.className='hex-plot-group';
    const grid=node('div');grid.className='hex-plot-grid';
    grid.append(...ids.map(id=>panel(id,true)));
    section.append(node('h2',title),grid);
    other.append(section);
  }
}
function renderHexPlot(host,{points,curve=[],line,xlabel},ylabel){
  const valid=points.filter(p=>Number.isFinite(p.x)&&Number.isFinite(p.y));
  if(!valid.length){const p=node('p','Недостаточно пригодных точек для этого графика.');p.className='plot-empty';host.append(p);return;}
  const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');
  svg.setAttribute('viewBox','0 0 680 330');svg.setAttribute('aria-label',ylabel+' — '+xlabel);svg.classList.add('plot');
  function draw(tag,attrs,label){
    const e=document.createElementNS(ns,tag);
    for(const [key,value] of Object.entries(attrs))e.setAttribute(key,value);
    if(label!==undefined)e.textContent=label;svg.append(e);return e;
  }
  const all=[...valid,...curve],xmin=Math.min(0,...all.map(p=>p.x));
  const xmax=Math.max(xmin+1e-14,...all.map(p=>p.x))*1.05;
  const lineY=line?[line.a*xmin+line.b,line.a*xmax+line.b]:[];
  const ymin=Math.min(0,...all.map(p=>p.y),...lineY),ymax=Math.max(0,...all.map(p=>p.y),...lineY);
  const spread=ymax-ymin||1,lo=ymin-spread*.04,hi=ymax+spread*.06;
  const x=v=>103+(v-xmin)*547/(xmax-xmin),y=v=>268-(v-lo)*243/(hi-lo);
  for(let j=0;j<=5;j++){
    const u=xmin+(xmax-xmin)*j/5,v=lo+(hi-lo)*j/5;
    draw('path',{class:'grid',d:'M'+x(u)+' 25V268 M103 '+y(v)+'H650'});
    draw('text',{x:x(u),y:289,'text-anchor':'middle'},plotNumber(u));
    draw('text',{x:97,y:y(v)+4,'text-anchor':'end'},plotNumber(v));
  }
  draw('path',{class:'axis',fill:'none',d:'M103 25V268H650'});
  if(line)draw('line',{class:'trend',x1:x(xmin),y1:y(lineY[0]),x2:x(xmax),y2:y(lineY[1])});
  if(curve.length)draw('polyline',{fill:'none',stroke:'#087c71','stroke-width':2,points:curve.map(p=>x(p.x)+','+y(p.y)).join(' ')});
  for(const p of valid){
    const label=[p.i===undefined?'':('Строка '+(p.i+1)),p.span,xlabel+' '+plotNumber(p.x),ylabel+' '+plotNumber(p.y),p.outside?'Экстраполяция':''].filter(Boolean).join(' · ');
    const attrs={class:'hex-plot-point',cx:x(p.x),cy:y(p.y),r:5,fill:p.outside?'#d96913':'#1c5395','aria-label':label};
    if(p.i!==undefined)Object.assign(attrs,{tabindex:0,role:'button','data-row':p.i});
    const dot=draw('circle',attrs);
    const tip=document.createElementNS(ns,'title');tip.textContent=label;dot.append(tip);
    if(p.i!==undefined){
      dot.onclick=()=>choose(p.i);
      dot.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose(p.i);}};
    }
  }
  draw('text',{class:'axis-label',x:365,y:319,'text-anchor':'middle'},xlabel);
  draw('text',{class:'axis-label',x:19,y:151,transform:'rotate(-90 19 151)','text-anchor':'middle'},ylabel);
  host.append(svg);
}
