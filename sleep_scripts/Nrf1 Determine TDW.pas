PROGRAM Determine TDW;
{classify wakefulness at theta dominated wakefulness according to 10.1073/pnas.1700983114}

USES
  strings;
TYPE
  epoch = RECORD
            state          : char;
            bin            : array[0..400] of single;
            EEGv,EMGv,temp : single;
          END;
VAR
  DIRf,TDWf                : text;
  SMOf                     : file of epoch;
  spec4                    : epoch;
  fln                      : string[3];
  str                      : string[2];
  t,x,m,strain,numM,lgth,
  wc                       : integer;
  fraction                 : single;
  vs                       : array[1..86500] of char;


PROCEDURE Calc_Fraction;
VAR
  theta,total,peakp       : single;
  hz,peakf                : integer;
BEGIN
  peakp:=0.0; FOR hz:=14 TO 60 DO IF spec4.bin[hz]>peakp THEN BEGIN peakp:=spec4.bin[hz]; peakf:=hz; END;  
  {determine frequency bin with max power between 3.5 and 15 Hz}
  IF (peakf>26) AND (peakf<48) THEN                                          {that max power has to be >6.5 Hz and <12 Hz. Previous larger range to avoid edge effects}
  BEGIN
    theta:=0.0; FOR hz:=(peakf-4) TO (peakf+4) DO theta:=theta+spec4.bin[hz];  {theta band 2.25 Hz broad adjusted according to peak frequency}
    total:=0.0; FOR hz:=14 TO 180 DO total:=total+spec4.bin[hz];               {total power 3.5 - 45 Hz}
    fraction:=theta/total;
  END
  ELSE fraction:=0.0;
END;


BEGIN
  assign(DIRf,'Nrf1_Directory.txt'); reset(DIRf);
  
  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);

    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5); fln:=trimright(fln5);
      assign(SMOf,fln+str+'.smo'); reset(SMOf);
      writeln(strain:3,m:4,fln+str:10);
	  
      t:=0;
      REPEAT
        inc(t);
        read(SMOf,spec4);
        vs[t]:=spec4.state;
        IF vs[t]='w' THEN {only in artefact-free waking}
        BEGIN
          Calc_Fraction;
          IF fraction>0.228 THEN vs[t]:='9'; {cut-off defined in 10.1073/pnas.1700983114}
        END
      UNTIL eof(SMOf);
      lgth:=t;
      close(SMOf);

      assign(TDWf,'D:\drp1\TDW2\'+fln+'.tdw'); rewrite(TDWf); 
	  {for each mouse a new 4day sequence of sleep-wake states is saved with '9' as TDW}
	  
      FOR t:=1 TO 3 DO writeln(TDWf,vs[t]);
      FOR t:=4 TO (lgth-3) DO
      BEGIN
        IF vs[t]='9' THEN
        BEGIN
          wc:=0;
          FOR x:=t-3 TO t+3 DO BEGIN IF (vs[x] in ['w','1','4']) THEN inc(wc); END;
          IF wc=6 THEN vs[t]:='w';
        END;
        IF (vs[t]='9') AND (vs[t-1] in ['n','2','5','r','3','6']) THEN vs[t]:='w';
        IF (vs[t]='9') AND (vs[t+1] in ['n','2','5','r','3','6']) THEN vs[t]:='w';
		{add context such that high theta in waking immediately before or after sleep  or 'isolated' 4s epochs of TDW are unlikely to be TDW; 
	     see 10.1073/pnas.1700983114}
        writeln(TDWf,vs[t]);
      END;
      FOR t:=lgth-2 TO lgth DO writeln(TDWf,vs[t]);
      close(TDWf);
    END;
  END;
  close(DIRf);
END.