PROGRAM LMA1h;
{locomotor activity (LMA), and waking and LMA/h wakefulness for 1h intervals over the 4-day recording}
{LMA recorded with a ClockLab at 1-min resolution independent from Somnologica}


USES
  strings;
VAR
  DIRf,ACTf,TDWf,OUTf      : text;
  fln5                     : string[5];
  fln                      : string;
  str                      : string[2];
  state                    : char;
  m,numM,h,t,strain,day,
  wak60                    : integer;
  local,act,act60,act1hw,
  time,wak60r              : single;
  

BEGIN
  assign(DIRf,'NRF Directory.txt'); reset(DIRf);
  assign(OUTf,'LMA 1h.out');  rewrite(OUTf);

  FOR strain:=1 TO 3 DO
  BEGIN
    readln(DIRf,str); readln(DIRf,numM);

    FOR m:=1 TO numM DO
    BEGIN
      readln(DIRf,fln5); fln:=trimright(fln5); writeln(str:3,m:3,fln:7);
      
	  assign(ACTf,fln+str+'.mldatx'); reset(ACTf); {open LMA files exported from ClockLab}
      assign(TDWf,fln+'.tdw'); reset(TDWf); 
	  {could also use SMO files to read sleep-wake states}
	  {but with the TDW files TDW and LMA could be directly compared, which is not done here}
     
      FOR h:=1 TO    3 DO readln(ACTf); {skip header}
      FOR h:=1 TO 3360 DO readln(ACTf); {skip 2 days + midnight to L-onset 2880+480 to snychronize with sleep-wake data}

      FOR day:=1 TO 4 DO
      BEGIN
        FOR h:=1 TO 24 DO
        BEGIN
		  act60:=0.0;
          FOR t:=1 TO 60 DO {read 1h of LMA}
          BEGIN
            readln(ACTf,local,local,local,act,local); {skip other variables in LMA file, keep activity}
            act60:=act60+act;
          END;

          wak60:=0;
          FOR t:=1 TO 900 DO {read 1h of waking}
          BEGIN
            readln(TDWf,state);
            IF state in ['w','1','4','9'] THEN inc(wak60);
          END;
          wak60r:=wak60; wak60r:=wak60r/900; {4sec to hours}
		  act1hw:=act60/wak60r; {LMA expressed per hour of wakefulness}
		  time:=h; time:=(time-0.5)+(day-1)*24; {for plotting data at interval midpoint over the 96h of the experiment}
		  writeln(OUTf,str:4,m:4,day:3,h:3,time:5:1,act60:12:3,wak60r*60:12:3,act1hw:12:3);
        END; {hour}
      END; {day}
	  
      close(ACTf); close(TDWf);
    END; {mouse}
  END; {strain}
  close(DIRf); close(OUTf);
END.