'use strict';

function workbookXmlText(value) {
  return String(value).replace(/[\u0000-\u0008\u000b\u000c\u000e-\u001f\ufffe\uffff]/g,'')
    .replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

function zipWorkbook(files) {
  const encoder=new TextEncoder();
  const crcTable=Array.from({length:256},(_,index)=>{
    let crc=index;
    for(let bit=0;bit<8;bit++)crc=(crc>>>1)^((crc&1)?0xedb88320:0);
    return crc;
  });
  const entries=[];
  const directory=[];
  let offset=0, directorySize=0;
  for(const [filename,xml] of Object.entries(files)){
    const name=encoder.encode(filename);
    const data=encoder.encode(xml);
    let crc=0xffffffff;
    for(const byte of data)crc=(crc>>>8)^crcTable[(crc^byte)&255];
    crc=(crc^0xffffffff)>>>0;

    const local=new Uint8Array(30+name.length);
    const header=new DataView(local.buffer);
    header.setUint32(0,0x04034b50,true);
    header.setUint16(4,20,true);
    header.setUint16(12,33,true);
    header.setUint32(14,crc,true);
    header.setUint32(18,data.length,true);
    header.setUint32(22,data.length,true);
    header.setUint16(26,name.length,true);
    local.set(name,30);
    entries.push(local,data);

    const central=new Uint8Array(46+name.length);
    const index=new DataView(central.buffer);
    index.setUint32(0,0x02014b50,true);
    index.setUint16(4,20,true);
    index.setUint16(6,20,true);
    index.setUint16(14,33,true);
    index.setUint32(16,crc,true);
    index.setUint32(20,data.length,true);
    index.setUint32(24,data.length,true);
    index.setUint16(28,name.length,true);
    index.setUint32(42,offset,true);
    central.set(name,46);
    directory.push(central);
    offset+=local.length+data.length;
    directorySize+=central.length;
  }
  const trailer=new Uint8Array(22);
  const footer=new DataView(trailer.buffer);
  footer.setUint32(0,0x06054b50,true);
  footer.setUint16(8,directory.length,true);
  footer.setUint16(10,directory.length,true);
  footer.setUint32(12,directorySize,true);
  footer.setUint32(16,offset,true);
  return new Blob([...entries,...directory,trailer],{
    type:'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
  });
}

function calibrationWorkbook(rows) {
  const xml='<?xml version="1.0" encoding="UTF-8"?>';
  const sheetNs='http://schemas.openxmlformats.org/spreadsheetml/2006/main';
  const officeRel='http://schemas.openxmlformats.org/officeDocument/2006/relationships';
  const packageRel='http://schemas.openxmlformats.org/package/2006/relationships';
  const headers=['Кислота','Концентрация, моль/л','Ожидаемая κ, мкСм/см','Текущая κ, мкСм/см'];
  const data=[headers,...rows].map((row,i)=>{
    const cells=row.map((value,j)=>{
      if(value===null)return '';
      const reference='ABCD'[j]+(i+1);
      if(typeof value==='string'){
        const style=i===0?' s="1"':'';
        return '<c r="'+reference+'" t="inlineStr"'+style+'><is><t xml:space="preserve">'+
          workbookXmlText(value)+'</t></is></c>';
      }
      const style=j>=2&&Math.abs(value)<100?' s="2"':'';
      return '<c r="'+reference+'"'+style+'><v>'+value+'</v></c>';
    }).join('');
    return '<row r="'+(i+1)+'">'+cells+'</row>';
  }).join('');
  return zipWorkbook({
    '[Content_Types].xml':xml+'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'+
      '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'+
      '<Default Extension="xml" ContentType="application/xml"/>'+
      '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'+
      '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'+
      '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/></Types>',
    '_rels/.rels':xml+'<Relationships xmlns="'+packageRel+'"><Relationship Id="rId1" Type="'+officeRel+
      '/officeDocument" Target="xl/workbook.xml"/></Relationships>',
    'xl/workbook.xml':xml+'<workbook xmlns="'+sheetNs+'" xmlns:r="'+officeRel+'"><bookViews><workbookView/></bookViews><sheets>'+
      '<sheet name="Ошибки калибровки" sheetId="1" r:id="rId1"/></sheets></workbook>',
    'xl/_rels/workbook.xml.rels':xml+'<Relationships xmlns="'+packageRel+'"><Relationship Id="rId1" Type="'+officeRel+
      '/worksheet" Target="worksheets/sheet1.xml"/><Relationship Id="rId2" Type="'+officeRel+
      '/styles" Target="styles.xml"/></Relationships>',
    'xl/styles.xml':xml+'<styleSheet xmlns="'+sheetNs+'"><numFmts count="1"><numFmt numFmtId="164" formatCode="0.0"/></numFmts>'+
      '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'+
      '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'+
      '<borders count="1"><border/></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'+
      '<cellXfs count="3"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'+
      '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'+
      '<xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/></cellXfs>'+
      '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>',
    'xl/worksheets/sheet1.xml':xml+'<worksheet xmlns="'+sheetNs+'"><dimension ref="A1:D'+(rows.length+1)+'"/>'+
      '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'+
      '<cols><col min="1" max="1" width="28" customWidth="1"/><col min="2" max="4" width="24" customWidth="1"/></cols>'+
      '<sheetData>'+data+'</sheetData><autoFilter ref="A1:D'+(rows.length+1)+'"/></worksheet>',
  });
}
