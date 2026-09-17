PROGRAM Delta_timecourse;

USES
  strings;
CONST
  day1=21600; day4=86400; {number of 4s epochs in 1 and in 4 days}
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
  tim,local,ref,del       : single;
  t,m,hz,strain,day,ld,
  start,intl,
  rest,h,nc,numM,numint,
  numi,s,e,so             : integer;
  vs                      : array[0..day4+1] of boolean;
  ref4                    : array[1..4] of single;


PROCEDURE readdel;
VAR
  tr : single;
BEGIN
  seek(EEGf,t-1); read(EEGf,spec4);
  local:=0.0; FOR hz:=4 TO 16 DO local:=local+spec4.bin[hz];
  local:=local*1000000/13;   {average power density in 1-4 Hz in microvolts square}
END;

PROCEDURE CALC_REF;
BEGIN
  FOR day:=1 TO 4 DO
  BEGIN
    nc:=0; ref4[day]:=0.0;
    t:=7199+(day-1)*day1;
    REPEAT
      inc(t);
      IF (vs[t]) AND (vs[t-1]) AND (vs[t+1]) THEN
      BEGIN
        inc(nc);
        readdel;
        ref4[day]:=ref4[day]+local;
      END;
    UNTIL t=10800+(day-1)*day1;
    ref4[day]:=ref4[day]/nc;
  END;
END;

PROCEDURE CALC_delta;
BEGIN
  nc:=0; del:=0.0; tim:=0.0;
  FOR h:=1 TO numint DO
  BEGIN
    IF h<=rest THEN numi:=intl+1 ELSE numi:=intl;
    REPEAT
      inc(t);
      IF (vs[t]) AND (vs[t-1]) AND (vs[t+1]) THEN
      BEGIN
        inc(nc);
        readdel;
        del:=del+local;
        tim:=tim+t;
      END;
    UNTIL (nc=numi) OR (t>=e);
    del:=del/nc; del:=del*100/ref; tim:=tim/nc; tim:=tim/900;
    writeln(OUTf,str:4,m:3,day:3,ld:3,h:4,tim:10:4,del:10:4,nc:8,ref:10:4);
    nc:=0; del:=0.0; tim:=0.0;
  END;
END;


BEGIN
  assign(DIRf,'Nrf1_Directory.txt'); reset(DIRf);
  assign(OUTf,'Nrf1_delta_timecourse.out'); rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);
    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5,so); fln:=trimright(fln5); writeln(str:3,m:3,fln:7);
      assign(EEGf,fln+str+'.smo'); reset(EEGf);

      {read 4 day sleep-wake state sequence}
      FOR t:=1 TO day4 DO
      BEGIN
        read(EEGf,spec4);
        IF spec4.state='n' THEN vs[t]:=true ELSE vs[t]:=false;
      END;
      vs[0]:=vs[1]; vs[day4+1]:=vs[day4];

      {calculate delta power in last 4h of the light periods to determine reference}
      CALC_REF;

      {calculate delta power time course over 4 days in percentiles}
      FOR day:=1 TO 4 DO
      BEGIN
        CASE day OF
          1,2 : ref:=(ref4[1]+ref4[2])/2;
          3   : ref:=(ref4[2]+ref4[4])/2;
          4   : ref:=ref4[4];
        END;

        FOR ld:=1 TO 2 DO
        BEGIN
          start:=(day-1)*day1;

          IF ld=1 THEN BEGIN numint:=12; s:=start;       END;   
		  IF ld=2 THEN BEGIN numint:=6;  s:=start+10800; END;
		  {in 12h light period 12 percentiles are used except after sleep deprivation when only 8 are used, the 12h dark periods are devided in 6 percentiles}
          e:=s+10800;
          IF (ld=1) AND (day=3) THEN BEGIN s:=so+start; numint:=8; END; {first percentile after sleep deprivation starts with sleep onset, read from the directory file}

          {determine number of 4s epochs per percentile}
          nc:=0; t:=s;
          REPEAT
            inc(t);
            IF (vs[t]) AND (vs[t-1]) AND (vs[t+1]) THEN inc(nc);   {only 4s epoch scored as artefact NREM sleep, i.e. 'n', and in between 2 other 'n' are considered}
          UNTIL t=e;
          intl:=nc DIV numint; rest:=nc MOD numint;

          t:=s;
          CALC_delta;
        END; {ld}
      END; {day}
      close(EEGf);
    END; {mouse}
  END; {strain}
  close(DIRf);
  close(OUTf);
END.