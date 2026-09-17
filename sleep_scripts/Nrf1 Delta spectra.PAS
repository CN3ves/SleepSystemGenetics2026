PROGRAM Delta Spectra;
{NREM sleep EEG spectra in the 0.5-4.5Hz range in the first recovery percentile after sleep deprivation}

USES
  strings;
CONST
  day1=21600; 
  day4=86400; {number of 4s per day and per experiment}
TYPE
  epoch = Record
            state          : char;
            bin            : array[0..400] of single;
            EEGv,EMGv,temp : single;
          end;
VAR
  DIRf,OUTf               : text;
  EEGf                    : file of epoch;
  spec4                   : epoch;
  str                     : string[2];
  fln5                    : string[5];
  fln                     : string;
  tim,ref,local,del,hzr   : single;
  localhz,delhz           : array[2..18] of single;
  t,m,hz,strain,day,ld,
  start,intl,
  rest,h,nc,numM,numint,
  numi,s,e,so             : integer;
  vs                      : array[0..day4+1] of boolean;
  ref4                    : array[1..4] of single;


PROCEDURE readdelhz;
BEGIN
  seek(EEGf,t-1); read(EEGf,spec4);
  FOR hz:=2 TO 18 DO localhz[hz]:=spec4.bin[hz]*1000000; {from V2/0.25Hz to microV2/0.25Hz}
END;

PROCEDURE CALC_REF;
BEGIN
  FOR day:=1 TO 4 DO
  BEGIN
    ref4[day]:=0.0; nc:=0;
    t:=7199+(day-1)*day1;
    REPEAT
      inc(t);
      IF (vs[t]) AND (vs[t-1]) AND (vs[t+1]) THEN
      BEGIN
        inc(nc);
        seek(EEGf,t-1); read(EEGf,spec4);
        local:=0; FOR hz:=4 TO 16 DO local:=local+spec4.bin[hz];
        local:=local*1000000/13;
        ref4[day]:=ref4[day]+local;
      END;
    UNTIL t=10800+(day-1)*day1;
    ref4[day]:=ref4[day]/nc;
  END;
END;

PROCEDURE CALC_delta;
BEGIN
  nc:=0; tim:=0.0; FOR hz:=2 TO 18 DO delhz[hz]:=0.0;
  h:=1; IF h<=rest THEN numi:=intl+1 ELSE numi:=intl;
  REPEAT
    inc(t);
    IF (vs[t]) AND (vs[t-1]) AND (vs[t+1]) THEN
    BEGIN
      inc(nc);
      readdelhz;
      FOR hz:=2 TO 18 DO delhz[hz]:=delhz[hz]+localhz[hz];
      tim:=tim+t;
    END;
  UNTIL (nc=numi) OR (t>=e);

  tim:=tim/nc; tim:=tim/900;
  FOR hz:=2 TO 18 DO
  BEGIN
    delhz[hz]:=delhz[hz]/nc;
    delhz[hz]:=delhz[hz]*100/ref;   {power density in all delta frequency bins are expressed to the same reference}
    hzr:=hz; hzr:=hzr/4;
    writeln(OUTf,str:4,m:3,tim:10:4,hzr:8:2,delhz[hz]:10:4,nc:8,ref:10:4);
  END;
END;


BEGIN
  assign(DIRf,'Nrf1_Directory.txt'); reset(DIRf);
  assign(OUTf,'Nrf1_delta_spec.out'); rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);

    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5,so); fln:=trimright(fln5); writeln(str:3,m:3,fln:7);
      assign(EEGf,fln+str+'.smo'); reset(EEGf);

      {read 4 day sleep-wake states}
      FOR t:=1 TO day4 DO
      BEGIN
        read(EEGf,spec4);
        IF spec4.state='n' THEN vs[t]:=true ELSE vs[t]:=false;
      END;
      vs[0]:=vs[1]; vs[day4+1]:=vs[day4];
      reset(EEGf);

      {calculate reference in the last 4h of the light periods, same as for used for delta power dynamics}
      CALC_REF;

      day:=3; ld:=1; ref:=(ref4[2]+ref4[4])/2;
      start:=(day-1)*day1;
      e:=start+10800; s:=so+start; numint:=8;

      {determine the number of 4s NREM sleep epochs for each of the 8 percentiles in the recovery day 1 light period}
      nc:=0; t:=s;
      REPEAT
        inc(t);
        IF (vs[t]) AND (vs[t-1]) AND (vs[t+1]) THEN inc(nc);
      UNTIL t=e;
      intl:=nc DIV numint; rest:=nc MOD numint;

      {calculate the delta spectra in the 1s percentile after sleep deprivation}
      t:=s;
      CALC_delta;

      close(EEGf);
    END; {mouse}
  END; {strain}
  close(DIRf);
  close(OUTf);
END.