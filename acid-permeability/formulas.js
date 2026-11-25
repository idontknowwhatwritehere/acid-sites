let pinnedFormula=null;
function formulaValue(value,unit=''){
  return Number.isFinite(value)?plotNumber(value)+(unit?' '+unit:''):'—';
}
function formulaTerm(symbol,label,value,unit,source,target=null){
  return {symbol,label,value,unit,source,target};
}
function renderFormulas(cal,run,row,isCal){
  const grid=el('formula-grid'),tooltip=el('formula-tooltip'),caption=el('formula-selection');
  tooltip.hidden=true;
  document.querySelectorAll('.formula-source').forEach(node=>node.classList.remove('formula-source','source-highlight'));
  const number=selected+1;
  const selectedRow=el('sheet').tBodies[0].rows[selected];
  const cell=index=>selectedRow?.cells[index];
  const param=id=>el('param-'+id);
  const calSource=state.sync?'Снимок i‑BuBn от '+new Date(state.sync.at).toLocaleString('ru-RU')+'.':'Исходные стандарты i‑BuBn.';
  const term=formulaTerm;
  const raw=term('κ','показание прибора',row?.raw,'мкСм/см','Строка '+number+' таблицы, столбец κ.',cell(2));
  const blank=term('κ₀','фон воды',cal.blank,'мкСм/см',calSource+' Настройка «Фон воды» для '+acid+'.');
  const signal=term('κ − κ₀','сигнал с поправкой',Number.isFinite(row?.raw)?row.raw-cal.blank:null,'мкСм/см','Из измеренной проводимости вычитается фон воды.');
  const a=term('a','наклон калибровки',cal.line?.a,cal.lambdas?'':'(мкСм/см)/(мкмоль/л)',calSource+' Регрессия по '+cal.points.length+' стандартам для '+acid+'.');
  const b=term('b','свободный член',cal.line?.b,'мкСм/см',calSource+' Пересечение калибровочной прямой с осью проводимости.');
  const x=term('X','ось калибровки',isCal?row?.x:row?.back?.x,cal.lambdas?'мкСм/см':'мкмоль/л',
    cal.lambdas?'Расчётная проводимость по сумме вкладов ионов.':'Расчётная концентрация H⁺ из равновесия диссоциации.');
  const fitted=term('ŷ','сигнал по прямой',cal.line&&Number.isFinite(x.value)?cal.line.a*x.value+cal.line.b:null,
    'мкСм/см','Сигнал с поправкой на фон, предсказанный калибровочной прямой.');
  const cards=[
    ['Поправка на фон',[signal,'=',raw,'−',blank]],
    ['Калибровочная прямая',[fitted,'=',a,'×',x,'+',b]]
  ];
  if(!isCal)cards.push(['Обратная калибровка',[x,'= (',signal,'−',b,') /',a]]);
  const h=term('[H⁺]','ионы водорода',isCal?row?.h:row?.back?.h,'мкмоль/л',
    cal.lambdas?'Найдено численно по проводимости и балансу зарядов.':'В опыте равно X из обратной калибровки.');
  const back=term('Cп','концентрация приёмника',row?.receiver,'М','Концентрация кислоты из обратной калибровки с учётом диссоциации.',cell(3));
  const c=isCal?term('C','концентрация стандарта',row?.c,'М',calSource+' Строка '+number+' таблицы стандартов.',cell(1)):back;
  const concentrationMicro=term(isCal?'Cµ':'Cп,µ','кислота в мкмоль/л',Number.isFinite(c.value)?c.value*1e6:null,'мкмоль/л','Концентрация в моль/л умножается на 10⁶.');
  const charge=term('q','средний заряд анионов',
    h.value>0&&cal.k?species(h.value,cal.k).charge:null,'',
    'q = α₁ + 2α₂ + 3α₃. Доли форм заданы константами диссоциации. K₁, K₂, K₃: '+(cal.k?.map(k=>formulaValue(k,'мкмоль/л')).join(', ')??'—')+'.');
  const oh=term('[OH⁻]','гидроксид-ионы',h.value>0?WATER.kw/h.value:null,'мкмоль/л','Автопротолиз воды: Kᵥ/[H⁺]. В этой шкале Kᵥ = 0,01 (мкмоль/л)².');
  cards.push(['Баланс зарядов',[h,'=',concentrationMicro,'×',charge,...(cal.water?['+',oh]:[])]]);
  if(cal.lambdas){
    const ions=h.value>0&&c.value>0?ionicConductivity(h.value,c.value,species(h.value,cal.k),cal.water,cal.lambdas):null;
    const conductance=term('Σλᵢcᵢ','сумма вкладов ионов',ions?(ions.value+ions.pure)*1000:null,'','Ионные проводимости λᵢ и рассчитанные концентрации cᵢ. λ° анионов: '+cal.lambdas.map(v=>fmt(v)).join(', ')+' См·см²/моль.');
    const conversion=term('1000','перевод единиц',1000,'','Множитель переводит произведение λᵢcᵢ в мкСм/см при cᵢ в мкмоль/л.');
    const pure=term('κводы','расчётный вклад воды',(WATER.h+WATER.oh)*Math.sqrt(WATER.kw)/1000,'мкСм/см','Проводимость чистой воды в модели автопротолиза.');
    cards.push(['Проводимость по ионам',[x,'=',conductance,'/',conversion,...(cal.water?['−',pure]:[])]]);
  }
  if(!isCal){
    const donor=term('Cд','концентрация донора',row?.donor,'М','Баланс вещества для выбранной строки.',cell(4));
    const vd=term('Vд','объём донора',run.vd,'мл','Входное поле «Объём донора».',param('vd'));
    const vr=term('Vп','объём приёмника',run.vr,'мл','Входное поле «Объём приёмника».',param('vr'));
    const cd0=term('Cд,0','донор в начале',run.cd,'М','Входное поле «C донора при t = 0».',param('cd'));
    const cr0=term('Cп,0','приёмник в начале',run.cr,'М','Входное поле «C приёмника при t = 0».',param('cr'));
    const area=term('A','площадь мембраны',run.area,'см²','Входное поле «Площадь мембраны».',param('area'));
    const delta=term('ΔC','градиент сейчас',row?.delta,'М','Концентрация донора минус концентрация приёмника.');
    const delta0=term('ΔC₀','начальный градиент',Number.isFinite(run.cd)&&Number.isFinite(run.cr)?run.cd-run.cr:null,'М','Разность введённых начальных концентраций.');
    const z=term('z','логарифм отношения',row?.z,'','Минус натуральный логарифм отношения текущего градиента к начальному.');
    const time=term('t','время замера',Number.isFinite(row?.t)?row.t*60:null,'с','Строка '+number+' таблицы, минуты умножены на 60.',cell(1));
    const beta=term('β','наклон z(t)',run.beta,'с⁻¹','МНК через начало координат по '+run.points.length+' точкам.');
    const sumTZ=term('Σ(tz)','сумма произведений',run.points.length?run.points.reduce((s,p)=>s+p.x*p.z,0):null,'с','Сумма t × z по точкам, включённым в регрессию.');
    const sumTT=term('Σ(t²)','сумма квадратов времени',run.points.length?run.sxx:null,'с²','Сумма t² по тем же строкам.');
    const p=term('P','общая проницаемость',run.result,'см/с','Наклон β, делённый на геометрический множитель двух камер.',el('result'));
    const n=term('n','перешедшее вещество',row?.transferred,'моль','Прирост количества кислоты в приёмнике.');
    const e=term('E','доля переноса',Number.isFinite(row?.transferred)&&run.vd>0&&run.cd>0?100*row.transferred/(run.vd/1000*run.cd):null,'%','Перешедшее количество относительно начального количества в доноре.');
    const interval=intervalDiagnostics(run).find(v=>v.i===selected);
    const flux=term('J','средний поток',interval?.flux,'моль/(см²·с)','Прирост n между соседними пригодными строками, делённый на площадь и длительность интервала.');
    const dn=term('Δn','прирост вещества',interval?.dn,'моль',interval?'Разность n в конце и начале интервала '+interval.span+'.':'Нужна предыдущая пригодная точка.');
    const dt=term('Δt','длительность интервала',interval?.dt,'с',interval?'Разность времени на интервале '+interval.span+'.':'Нужна предыдущая пригодная точка.');
    const pi=term('Pинт','проницаемость интервала',interval?.p,'см/с','Оценка по изменению z между соседними пригодными замерами.');
    const dz=term('Δz','изменение z',interval?.dz,'','z в конце интервала минус z в начале.');
    const geo=term('G','геометрический множитель',run.factor,'см⁻¹','A × (1/Vд + 1/Vп) при объёмах в см³.');
    const err=row?.eligible&&run.beta!==null?row.z-run.beta*row.t*60:null;
    const residual=term('ε','остаток',err,'','Фактическое значение z минус прогноз βt.');
    const leverage=row?.eligible&&run.sxx>0?(60*row.t)**2/run.sxx:null;
    const mse=run.sse!==null&&run.points.length>1?run.sse/(run.points.length-1):null;
    const cook=err!==null&&leverage!==null&&leverage<1-1e-12&&mse>0?err**2*leverage/(mse*(1-leverage)**2):null;
    const hInfluence=term('hᵢ','рычаг точки',leverage,'','tᵢ² / Σ(t²) для регрессии через начало.');
    const mseTerm=term('MSE','дисперсия остатков',mse,'','SSE / (n − 1), где n — число включённых точек и оценивается один наклон.');
    const sse=term('SSE(P)','сумма квадратов',run.result!==null?run.sse:null,'','Сумма по всем включённым строкам. Показано значение при оценённом P.');
    const startP=hexPlotData(run,'start').points.find(v=>v.i===selected)?.y;
    cards.push(['Баланс вещества в доноре',[donor,'=',cd0,'− (',vr,'/',vd,') × (',back,'−',cr0,')']],
      ['Перенесённое количество',[n,'= (',vr,'/ 1000) × (',back,'−',cr0,')']],
      ['Степень переноса',[e,'= 100 ×',n,'/ ((',vd,'/ 1000) ×',cd0,')']],
      ['Разность концентраций',[delta,'=',donor,'−',back]],
      ['Логарифмическая линеаризация',[z,'= −ln(',delta,'/',delta0,')']],
      ['Наклон общей регрессии',[beta,'=',sumTZ,'/',sumTT]],
      ['Остаток регрессии',[residual,'=',z,'−',beta,'×',time]],
      ['Проницаемость серии',[p,'=',beta,'/ [',area,'× (1/',vd,'+ 1/',vr,')]']],
      ['Интервальная проницаемость',[pi,'=',dz,'/ (',geo,'×',dt,')']],
      ['Средний поток',[flux,'=',dn,'/ (',area,'×',dt,')']],
      ['Сумма квадратов',[sse,'= Σ(',z,'−',p,'×',geo,'×',time,')²']],
      ['Расстояние Кука',[term('Dᵢ','влияние точки',cook,'','Регрессия через начало координат с одним параметром.'),'= (',residual,'² /',mseTerm,') ×',hInfluence,'/ (1 −',hInfluence,')²']],
      ['Сдвиг начала участка',[term('P(tнач)','проницаемость участка',startP,'см/с','Наклон по последующим точкам после вычитания t и z выбранного начала.'),'=',
        term('βуч','наклон участка',Number.isFinite(startP)?startP*run.factor:null,'с⁻¹','Σ[(t−tнач)(z−zнач)] / Σ[(t−tнач)²].'),'/',geo]]);
    renderSteps(number,[
      ['Измерение и поправка на фон',[raw,blank,signal]],
      ['Обратная калибровка',[a,b,x,back]],
      ['Баланс вещества',[cr0,vr,n,vd,cd0,donor,e]],
      ['Градиент и логарифм',[delta0,delta,z]],
      ['Интервал до выбранной точки',[dt,dn,flux,pi]],
      ['Регрессия по серии',[sumTZ,sumTT,beta,geo,p,sse,residual]]
    ]);
  } else renderSteps(number,[['Исходный стандарт',[c,raw,blank]],['Расчёт оси X',[h,charge,x]],['Калибровочная прямая',[signal,a,b]]]);
  const previous=pinnedFormula;
  let matched=false;
  const entries=cards.map(([title,parts])=>{
    const card=node('article');card.className='formula-card';
    const strip=node('div');strip.className='formula-terms';
    for(const [i,part] of parts.entries()){
      if(typeof part==='string'){strip.append(node('span',part));continue;}
      const id=title+'-'+i;
      const button=node('button');button.className='formula-term';button.type='button';
      button.setAttribute('aria-pressed',String(id===previous));
      button.append(node('span',part.symbol),node('small',part.label));
      const description=()=>part.symbol+' · '+part.label+' · '+formulaValue(part.value,part.unit)+'. '+part.source;
      if(id===previous){caption.textContent=description();matched=true;}
      const hide=()=>{tooltip.hidden=true;button.removeAttribute('aria-describedby');part.target?.classList.remove('source-highlight','formula-source');};
      const show=()=>{
        tooltip.textContent=description();tooltip.hidden=false;button.setAttribute('aria-describedby','formula-tooltip');
        part.target?.classList.add('source-highlight','formula-source');
        const rect=button.getBoundingClientRect();
        const left=Math.max(8,Math.min(rect.left,innerWidth-tooltip.offsetWidth-8));
        const above=rect.top-tooltip.offsetHeight-6;
        tooltip.style.left=left+'px';
        tooltip.style.top=(above>=8?above:Math.max(8,Math.min(rect.bottom+6,innerHeight-tooltip.offsetHeight-8)))+'px';
      };
      button.onpointerenter=show;button.onpointerleave=hide;
      button.onfocus=show;button.onblur=hide;
      button.onkeydown=event=>{if(event.key==='Escape')hide();};
      button.onclick=()=>{
        pinnedFormula=id;
        grid.querySelector('[aria-pressed=true]')?.setAttribute('aria-pressed','false');
        button.setAttribute('aria-pressed','true');caption.textContent=description();show();
      };
      strip.append(button);
    }
    card.append(node('h3',title),strip);return card;
  });
  if(!matched)caption.textContent='Выбери член формулы, чтобы закрепить значение и источник.';
  grid.replaceChildren(...entries);
}
function renderSteps(number,steps){
  el('step-row').textContent='Строка '+number;
  el('calculation-steps').replaceChildren(...steps.map(([label,terms])=>{
    const li=node('li');
    li.append(node('strong',label),node('p',terms.map(t=>t.symbol+' = '+formulaValue(t.value,t.unit)).join(' · ')));
    return li;
  }));
}
