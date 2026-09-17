PROGRAM TDW timecourse;
{calculate hourly values for time spent in wkaing, TDW, and %TDW of waking}

USES
  strings;
VAR
  DIRf,OUTf,TDWf       : text;
  fln5                 : string[5];
  str                  : string[2];
  vs                   : char;
  m,numM,h,t,wc,tc,ac,
  strain,day,start     : integer;
  wr,tr,ar,time,cor    : single;


BEGIN
  assign(DIRf,'Nrf1 directory.txt'); reset(DIRf);
  assign(OUTf,'Nrf1 TDW1h.out'); rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);

    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5); fln:=trimright(fln5);
      assign(TDWf,'fln+'.tdw'); reset(TDWf); {files made with Nrf1 determine TDW.pas}

      FOR day:=1 TO 4 DO
      BEGIN
        FOR h:=1 TO 24 DO
        BEGIN
          wc:=0; tc:=0; ac:=0;
          FOR t:=1 TO 900 DO
          BEGIN
            readln(TDWf,vs);
            CASE vs OF
              'w'     : inc(wc);
              '9'     : inc(tc);
              '1','4' : inc(ac);
            END;
          END;

          wr:=wc; tr:=tc; ar:=ac;
          IF wr>0.0 THEN
          BEGIN
            cor:=tr/(wr+tr)*ar; 
			{as TDW is only determined in artefact-free waking,
			the hourly time spent in TDW is compensated proportionally to the number of 4s of waking with artefacts within a given hour}
            tr:=tr+cor;
            ar:=ar-cor;
            wr:=wr+tr+ar;
            ar:=tr*100/wr; {TDW as % of waking, sorry to reuse the artefact counter}
            wr:=wr/15; {waking min/h}
            tr:=tr/15; {TDW min/h}
          END;
          IF ar>100.0 THEN BEGIN write(wc:5,tc:5,ac:5); readln; END;
          time:=h+(day-1)*24;
          time:=time-0.5; {for plotting at interval midpoint}
          writeln(OUTf,str:3,m:4,day:3,h:3,time:6:1,wr:10:4,tr:10:4,ar:10:4);
        END; {h}
      END; {day}
      close(TDWf);
    END;
  END;
  close(DIRf);
  close(OUTf);
END.