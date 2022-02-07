#pragma rtGlobals=1		// Use modern global access method.
#include "ethographer"
#include<power spectral density>
#include<GIS Utilities>
#include <RosePlot>
#include <Peak AutoFind>
#include <Multi-peak fitting 2.0>
#include <FilterDialog> menus=0
#include <All IP Procedures>
#include <Image Saver>

//////////////////////////////////////////////////////////////////////
//////////     すべてをまとめた関数　ALL    //////////
//////////////////////////////////////////////////////////////////////
function ALL()
variable p=13
string ID, spname="ss", utmlat, utmlon, spd, direc, turn, startp, endp, rad
do
	ID = num2istr(p)
	wave lat = $(spname+ID+"lat"), lon = $(spname+ID+"lon")
	utmlat=spname+ID+"utmlat"; utmlon=spname+ID+"utmlon"; spd=spname+ID+"spd"; direc=spname+ID+"dir"
	turn=spname+ID+"turn"; startp=spname+ID+"start"; endp=spname+ID+"end"; rad=spname+ID+"rad"
	if(waveexists($(spname+ID+"spd"))==0)
		flightmask_func(lat, lon, spname, ID, utmlat, utmlon, spd, direc, turn, startp, endp, rad)
	endif
	wave flightstart = $(spname+ID+"start"), flightend = $(spname+ID+"end")
	if(waveexists($(spname+ID+"sp"))==0)
		spep(flightstart, flightend, spname, ID)
	endif
	wave radwave = $(spname+ID+"rad"), startwave = $(spname+ID+"start"), endwave = $(spname+ID+"end"), spdwave = $(spname+ID+"spd")
//	if(waveexists($(spname+ID+"sp"))==0)
//		commuting_func(radwave, startwave, endwave, 0, 0, 0, 0, spdwave, spname, ID)
//	endif
	wave sp = $(spname+ID+"sp"), ep = $(spname+ID+"ep"), radwave = $(spname+ID+"rad")
	if(waveexists($(spname+ID+"windspd"))==0)
		windestimate(sp, ep, spname, ID, lat, lon, spdwave, radwave)//////////choose graph number//////////choose air speed
	endif
	if(waveexists($(spname+ID+"windspd_ca"))==0)
		windestimate_const_airspeed(sp, ep, spname, ID, lat, lon, spdwave, radwave)//////////choose graph number//////////choose air speed
	endif
//	wave radianwave = $(spname+ID+"windrad")
//	if(waveexists($(spname+ID+"bwang"))==0)
//		pathwindangle(radwave, spdwave, sp, ep, radianwave, spname, ID)
//	endif
//	wave utmlatwave=$(spname+ID+"utmlat"), utmlonwave=$(spname+ID+"utmlon"), windspd=$(spname+ID+"windspd")
//	if(waveexists($(spname+ID+"st"))==0)
//		st(utmlatwave, utmlonwave, sp, ep, spdwave, windspd, spname, ID) 
//	endif
//	wave utmlatwave=$(spname+ID+"utmlat"), utmlonwave=$(spname+ID+"utmlon")
//	if(waveexists($(spname+ID+"smthst"))==0)
//		smooth_func(utmlatwave, utmlonwave, sp, ep, spdwave, spname, ID, windspd)
//	endif
p+=1
while(p<=13)
end

function ALLswctwc()
variable p
string ID, spname="wan", windspd, bwang
p=1
do
ID=num2str(p)
wave windspdwave=$(spname+ID+"windspd"), bwangwave=$(spname+ID+"bwang")
swctwc(windspdwave,bwangwave,spname,ID)
p+=1
while(p<=4)
end

function compile()
//wave ss1bwang, ss2bwang, ss3bwang, ss4bwang, ss5bwang, ss6bwang, ss7bwang
//wave ss11bwang, ss12bwang, ss13bwang, ss14bwang, ss15bwang, ss16bwang
//wave ss31bwang, ss32bwang, ss33bwang
wave wan1bwang, wan2bwang, wan3bwang, wan4bwang
//concatenate/O {ss1bwang, ss2bwang, ss3bwang, ss4bwang, ss5bwang, ss6bwang, ss7bwang}, allbwang1
//concatenate/O {ss11bwang, ss12bwang, ss13bwang, ss14bwang, ss15bwang, ss16bwang}, allbwang2
//concatenate/O {ss31bwang, ss32bwang, ss33bwang}, allbwang3
//concatenate/O {allbwang1, allbwang2, allbwang3}, ALLbwang
concatenate/O {wan1bwang, wan2bwang, wan3bwang, wan4bwang}, WANbwang
end
//////////////////////////////////////////////////////////////////////
//////////     飛行区間を抽出して開始点と終了点を導く　flightmask    //////////
//////////////////////////////////////////////////////////////////////
function flightmask_func(ss1lat,ss1lon, spname, num, utmlat, utmlon, spd, direc, turn, startp, endp, rad)   /////   出力wave utmlat, utmlon, spd, dir, turn, rad, start, end, flight
wave ss1lat, ss1lon
string spname, num, utmlat, utmlon, spd, direc, turn, startp, endp, rad
	transform2utm(ss1lat, ss1lon)
calcspddirturn(utm_y, utm_x)
rename utm_y $utmlat; rename utm_x $utmlon
rename utm_x_spd $spd; rename utm_x_dir $direc; rename utm_x_turn $turn
wave spdwave = $(spname+num+"spd")
spdwave = spdwave * 1000 / 3600
makemask(spdwave, low=4, high=inf, mname="flight")
wave flight
duplicate flight flightrev
maskinverse(flightrev)
masklength(flightrev, low=0, high=5)
maskcalculation("flight", "flightrev", "OR", spname+num+"flight")
wave flightwave = $(spname+num+"flight")
masklength(flightwave, low=600, high=inf)
killwaves flight flightrev
insertpoints numpnts(flightwave)-1,1,flightwave//終了点からのマスクを判定するために仮に0を入れる
insertpoints 0,1,flightwave//開始点からのマスクを判定するために仮に0を入れる
statmask(ss1lat, flightwave)
getstatdata(stat_result, "stat", "startrow")
getstatdata(stat_result, "stat", "endrow")
wave startrow_stat, endrow_stat//仮に入れた0の分を引く
startrow_stat = startrow_stat-1; endrow_stat=endrow_stat-1
rename startrow_stat $startp; rename endrow_stat $endp
killwaves stat_result
wave direction = $(spname+num+"dir")
duplicate direction $rad
wave radian = $(spname+num+"rad")
radian = direction * pi/180
//wave ss1start, ss1end
//edit ss1start ss1end
wave turnwave = $(spname+num+"turn")
killwaves turnwave
end

function flightmask()

flightmask_func(ss1lat, ss1lon, "ss", "1", "ss1utmlat", "ss1utmlon", "ss1spd", "ss1dir", "ss1turn", "ss1start", "ss1end", "ss1rad")

end

////////////////////////////////////////////////////////////////////////////////////////////
//////////          　進行方向ベクトルの平均からcommuting flight だけを抽出する関数　commuting()         /////////
////////////////////////////////////////////////////////////////////////////////////////////
function commuting_func(ss1rad, ss1start, ss1end, p, q, r, n, ss1spd, spname, num)//ss1vecシリーズ(time, length, start, end)と、ss1commuteができる
wave ss1rad, ss1start, ss1end, ss1spd
//ss1vectimeはあらかじめ大きい値を1つ入れたウェーブを作っておく。時刻ウェーブを正しく表示するため。
variable p,q,r,n
string spname, num
make/O/N=10000 $(spname+num+"vecstart"), $(spname+num+"vecend"), $(spname+num+"veclength")
wave vecstart = $(spname+num+"vecstart"), vecend = $(spname+num+"vecend"), veclength = $(spname+num+"veclength")
veclength=nan //進行方向のベクトル平均データを格納するウェーブ
	r=0
	for(q=0;q<numpnts(ss1start);q+=1)
	for(p=0;p<(ss1end[q] - ss1start[q] - 300 - 60)/30;p+=1)//飛行開始後と終了前の60秒は含めない
	duplicate/R=[ss1start[q]+60+p*30, ss1start[q]+60+p*30+300] ss1rad ss1rad_ref //5分の窓を30秒秒間隔で移動させて進行方向ベクトルの長さ平均を計算する
	duplicate ss1rad_ref ss11; ss11 = 1 //ベクトルの長さは1で統一
	concatenate {ss1rad_ref, ss11}, ss1rad_vec
	statscircularmeans/Q ss1rad_vec
		wave w_circularmeans
		veclength[r] = w_circularmeans[1]
		//vectime[r] = ss1time[ss1start[q]+p*30+150]
		vecstart[r] = ss1start[q]+60+p*30; vecend[r] = ss1start[q]+60+p*30+300//上に合わせて区間の開始と終了を決める
		killwaves ss1rad_ref, ss11, ss1rad_vec, w_circularmeans
	r+=1
	endfor
	endfor
		p=0//nanを消す
		do
		if(numtype(veclength[p])==2)
		deletepoints p,1,veclength, vecstart, vecend
		else
		p+=1
		endif
		while(p<numpnts(veclength))
//ss1commuteを作る。5分間の進行方向ベクトルの長さ平均が0.7を超える部分のみを抽出。
duplicate/O ss1spd, $(spname+num+"commute")
wave commute = $(spname+num+"commute"); commute = 0
	for(p=0; p<numpnts(veclength); p+=1)
	if(veclength[p]>=0.70)
	commute[vecstart[p], vecend[p]] = 1
	endif
	endfor
	MakeMask(commute, Low=0.5, High=inf, mName=spname+num+"commutemask")
	killwaves commute
wave commutemask = $(spname+num+"commutemask")//comspとcomepを作る
	StatMask(commutemask, commutemask)
	GetStatData(Stat_Result, "Stat", "startRow")
	GetStatData(Stat_Result, "Stat", "endRow")
	wave startrow_stat, endrow_stat
		make/O/N=3000 $(spname+num+"sp"), $(spname+num+"ep")
		wave comsp = $(spname+num+"sp"), comep = $(spname+num+"ep")
		comsp = nan; comep = nan
		p = 0; n = 0
		do
		q = 0
		do
		comsp[p] = startrow_stat[n] + 300*q//5分の窓を5分ずつずらしていく
		comep[p] = startrow_stat[n] + 300*q + 300
		p+=1; q+=1
		while(q<((endrow_stat[n]-startrow_stat[n])-300)/300)
		n+=1
		while(n<numpnts(startrow_stat))
			p=0
			do
			if(numtype(comsp[p])==2)
			deletepoints p, 1, comsp, comep
			else
			p+=1
			endif
			while(p<numpnts(comsp))
	killwaves stat_result, startrow_stat, endrow_stat
killwaves vecstart, vecend, veclength
end

function commuting()
variable p=1
string ID, spname="ss"
do
ID = num2istr(p)
wave rad = $(spname+ID+"rad"), startp = $(spname+ID+"start"), endp = $(spname+ID+"end"), spd = $(spname+ID+"spd")
commuting_func(rad, startp, endp, 0, 0, 0, 0, spd, "ss", ID)
p+=1
while(p<=2)
end

/////////////////////////////////////////////////////////////////////////////
//////////          5分間隔の開始点と終了点を比較する spep()         //////////
/////////////////////////////////////////////////////////////////////////////
function spep(ss1start, ss1end, spname, num)
wave ss1start, ss1end
string spname, num
variable p, q, r
make/N=3000/O $(spname+num+"sp"), $(spname+num+"ep")
wave sp=$(spname+num+"sp"), ep=$(spname+num+"ep")
sp=nan; ep=nan
p=0; q=0
do
r=0
do
	sp[p] = ss1start[q] + 60 + r*300
	ep[p] = ss1start[q] + 60 + r*300 + 300
	p+=1; r+=1
while(r<(ss1end[q]-ss1start[q]-120-300)/300)
	q+=1
while(q<numpnts(ss1start))
p=0
do
	if(numtype(sp[p])==2)
	deletepoints p,1,sp,ep
	else
	p+=1
	endif
while(p<numpnts(sp))
end



/////////////////////////////////////////////////////////////////////////////
//////////          移動飛行している区間の風を5分間隔で推定する windestimate()         //////////
/////////////////////////////////////////////////////////////////////////////
function windestimate(ss1sp, ss1ep, spname, num, ss1lat, ss1lon, ss1spd, ss1rad)   /////   出力
wave ss1sp, ss1ep, ss1lat, ss1lon, ss1spd, ss1rad
string spname, num
variable p
make/N=5000/O/D $(spname+num+"windrad"), $(spname+num+"windspd"), $(spname+num+"lonave")
make/N=5000/O/D $(spname+num+"latave"), $(spname+num+"sinaic"), $(spname+num+"linaic"), $(spname+num+"sinlin")
wave windrad=$(spname+num+"windrad"), windspd=$(spname+num+"windspd"), lonave=$(spname+num+"lonave")
wave latave=$(spname+num+"latave"), sinaic=$(spname+num+"sinaic"), linaic=$(spname+num+"linaic"), sinlin=$(spname+num+"sinlin")
windrad=nan; windspd=nan; lonave=nan; latave=nan; sinaic=nan; linaic=nan; sinlin=nan
p=0
do
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1lat latcut
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1lon loncut
	lonave[p] = mean(loncut); latave[p] = mean(latcut)
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1spd speed
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1rad angle
		display speed vs angle; setaxis bottom 0, 2*pi
		Make/D/N=3/O W_coef; W_coef[0] = {10,0,0}//////////     choose air speed / or choose not to constrain air speed
		FuncFit/Q=1/L=3600/X=1/NTHR=0 sinusoidal W_coef  speed /X=angle /D /R //////////     sine fitting
		WaveStats/Q fit_speed
		windrad[p] = pi/2 - V_maxloc
		windspd[p] = (V_max - V_min)/2
		wave Res_speed
		duplicate/O Res_speed res; res = Res_speed ^ 2
		sinaic[p] = numpnts(res)*((log(2*pi*sum(res)/numpnts(res))/log(e))+1)+2*(1+3)
		dowindow/K graph1 //////////     choose graph number
	killwaves W_coef W_sigma fit_speed res_speed res
		display speed vs angle; setaxis bottom 0, 2*pi
		K1 = 0;
		CurveFit/Q=1/L=360 /X=1/H="01"/NTHR=0 line  speed /X=angle /D /R //////////     line fitting
		wave Res_speed
		duplicate/O Res_speed res
		res = Res_speed ^ 2
		linaic[p] = numpnts(res)*((log(2*pi*sum(res)/numpnts(res))/log(e))+1)+2*(0+2)
	if(sinaic[p]>=linaic[p]-2) //////////     compare sine and line
	sinlin[p] = 1
	else
	sinlin[p] = 0
	endif
		dowindow/K graph1 //////////     choose graph number
	if(windspd[p] > 100 || windspd[p] < 0.01)
	windrad[p]=0; windspd[p]=0; lonave[p]=0; latave[p]=0
	sinaic[p]=0; linaic[p]=0; sinlin[p]=1
	endif
killwaves speed angle latcut loncut W_coef fit_speed W_sigma res_speed res
p+=1
while (p < numpnts(ss1sp))
edit windrad windspd latave lonave ss1sp ss1ep sinaic linaic sinlin
variable s
	s=0
	do
	if(numtype(windrad[S])==2)
	deletepoints S, 1, windrad; deletepoints S, 1, windspd; deletepoints S, 1, latave
	deletepoints S, 1, lonave; deletepoints S, 1, sinaic; deletepoints S, 1, linaic; deletepoints S, 1, sinlin
	else
	S+=1
	endif
	while(S<numpnts(windrad))
	S=0
	do
	if(sinlin[S]==1)
	windrad[s]=nan; windspd[s]=nan; latave[s]=nan; lonave[s]=nan; sinaic[s]=nan; linaic[s]=nan; sinlin[s]=nan;
	endif
	s+=1
	while(S<numpnts(sinlin))
//	killwaves sinlin, sinaic, linaic
//////////　　　windvectorデータを作成する
duplicate/O windspd a; a = windspd * 10
concatenate {a, windrad}, $(spname+num+"windvec")
killwaves a
End

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          移動飛行している区間の風を5分間隔で推定する（対気速度一定） windestimate_const_airspeed()         //////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
function windestimate_const_airspeed(ss1sp, ss1ep, spname, num, ss1lat, ss1lon, ss1spd, ss1rad)   /////   出力
wave ss1sp, ss1ep, ss1lat, ss1lon, ss1spd, ss1rad
string spname, num
variable p
make/N=5000/O/D $(spname+num+"windrad_ca"), $(spname+num+"windspd_ca")
make/N=5000/O/D $(spname+num+"sinaic_ca"), $(spname+num+"linaic_ca"), $(spname+num+"sinlin_ca")
wave windrad=$(spname+num+"windrad_ca"), windspd=$(spname+num+"windspd_ca")
wave sinaic=$(spname+num+"sinaic_ca"), linaic=$(spname+num+"linaic_ca"), sinlin=$(spname+num+"sinlin_ca")
windrad=nan; windspd=nan; sinaic=nan; linaic=nan; sinlin=nan
p=0
do
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1lat latcut
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1lon loncut
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1spd speed
	duplicate/O/R = [ss1sp[p],ss1ep[p]] ss1rad angle
		display speed vs angle; setaxis bottom 0, 2*pi
		Make/D/N=3/O W_coef; W_coef[0] = {10,0,0}//////////     choose air speed / or choose not to constrain air speed
		FuncFit/Q=1/L=3600/X=1/H="100"/NTHR=0 sinusoidal W_coef  speed /X=angle /D /R //////////     sine fitting
		WaveStats/Q fit_speed
		windrad[p] = pi/2 - V_maxloc
		windspd[p] = (V_max - V_min)/2
		wave Res_speed
		duplicate/O Res_speed res; res = Res_speed ^ 2
		sinaic[p] = numpnts(res)*((log(2*pi*sum(res)/numpnts(res))/log(e))+1)+2*(1+3)
		dowindow/K graph1 //////////     choose graph number
	killwaves W_coef W_sigma fit_speed res_speed res
		display speed vs angle; setaxis bottom 0, 2*pi
		K1 = 0;
		CurveFit/Q=1/L=360 /X=1/H="01"/NTHR=0 line  speed /X=angle /D /R //////////     line fitting
		wave Res_speed
		duplicate/O Res_speed res
		res = Res_speed ^ 2
		linaic[p] = numpnts(res)*((log(2*pi*sum(res)/numpnts(res))/log(e))+1)+2*(0+2)
	if(sinaic[p]>=linaic[p]-2) //////////     compare sine and line
	sinlin[p] = 1
	else
	sinlin[p] = 0
	endif
		dowindow/K graph1 //////////     choose graph number
	if(windspd[p] > 100 || windspd[p] < 0.01)
	windrad[p]=0; windspd[p]=0
	sinaic[p]=0; linaic[p]=0; sinlin[p]=1
	endif
killwaves speed angle latcut loncut W_coef fit_speed W_sigma res_speed res
p+=1
while (p < numpnts(ss1sp))
edit windrad windspd ss1sp ss1ep// sinaic linaic sinlin
variable s
	s=0
	do
	if(numtype(windrad[S])==2)
	deletepoints S, 1, windrad; deletepoints S, 1, windspd
	deletepoints S, 1, sinaic; deletepoints S, 1, linaic; deletepoints S, 1, sinlin
	else
	S+=1
	endif
	while(S<numpnts(windrad))
	S=0
	do
	if(sinlin[S]==1)
	windrad[s]=nan; windspd[s]=nan; sinaic[s]=nan; linaic[s]=nan; sinlin[s]=nan;
	endif
	s+=1
	while(S<numpnts(sinlin))
	killwaves sinlin, sinaic, linaic
//////////　　　windvectorデータを作成する
duplicate/O windspd a; a = windspd * 10
concatenate {a, windrad}, $(spname+num+"windvec_ca")
killwaves a
End


////////////////////////////////////////////////////////////////////////////////////////////////
//////////          それぞれの区間の鳥の進行方向と風とのなす角度を計算する関数　pathwindangle()         //////////
////////////////////////////////////////////////////////////////////////////////////////////////
// 進行方向ベクトルの平均方向を鳥の進行方向とみなした場合
function pathwindangle(ss1rad, ss1spd, ss1sp, ss1ep, ss1windrad, spname, num)   /////   出力wave vechead, vecwindang
wave ss1rad, ss1spd, ss1sp, ss1ep, ss1windrad
string spname, num
variable p
make/O/N=(numpnts(ss1sp)) $(spname+num+"head"), $(spname+num+"length")
wave head=$(spname+num+"head"), length=$(spname+num+"length")
p=0
do
	duplicate/O/R=[ss1sp[p], ss1ep[p]] ss1rad radcut
	duplicate/O/R=[ss1sp[p], ss1ep[p]] ss1spd spdcut
	concatenate/O {radcut, spdcut}, veccut
	statscircularmeans/Q veccut
	wave w_circularmeans
	head[p] = w_circularmeans[2]
	length[p] = w_circularmeans[1]
	if(head[p]>=0 && head[p] <pi)//ベクトルの向きを風ベクトル（x軸が0度）と合わせる
	head[p] = pi/2 - head[p]
	else
	head[p] = (5*pi)/2 - head[p]
	endif 
	p+=1
while(p<numpnts(ss1sp))
killwaves length
variable q
duplicate/O ss1windrad ss1r; ss1r = ss1windrad + 3/2*pi; duplicate/O head ss1h; ss1h = head + 3/2*pi
duplicate/O ss1r $(spname+num+"windang")
wave windang=$(spname+num+"windang")
q=0
do
	if(ss1r[q]-ss1h[q] >= 0 && ss1r[q]-ss1h[q] <=pi)
	windang[q] = ss1r[q]-ss1h[q]
	endif
	if(ss1r[q]-ss1h[q] < 0 && ss1r[q]-ss1h[q] >=-pi)
	windang[q] = ss1h[q]-ss1r[q]
	endif
	if(ss1r[q]-ss1h[q] > pi && ss1r[q]-ss1h[q] <=2*pi)
	windang[q] = 2*pi-(ss1r[q]-ss1h[q])
	endif
	if(ss1r[q]-ss1h[q] < -pi && ss1r[q]-ss1h[q] >=-2*pi)
	windang[q] = 2*pi+(ss1r[q]-ss1h[q])
	endif
	q+=1
	while(q<numpnts(ss1r))
killwaves ss1r ss1h
duplicate/O windang $(spname+num+"bwang")
wave bwang=$(spname+num+"bwang")
for(p=0; p<numpnts(windang); p+=1)
if(windang[p] > pi)
bwang[p] = 2*pi - windang[p]
endif
endfor
killwaves windang
killwaves radcut spdcut veccut w_circularmeans
bwang = bwang *180/pi
end

/////////////////////////////////////////////////////////////////////////////////////////////
//////////          5分区間の開始点と終了点を結んだ直線と実経路の長さの割合を計算する関数　st         //////////
/////////////////////////////////////////////////////////////////////////////////////////////
function st(ss1utmlat, ss1utmlon, ss1sp, ss1ep, ss1spd, ss1windspd, spname, num)     /////   出力wave st, stdis
wave ss1utmlat, ss1utmlon, ss1sp, ss1ep, ss1spd, ss1windspd
string spname, num
variable p
duplicate/O ss1sp $(spname+num+"stdis"); duplicate/O ss1sp $(spname+num+"st"); duplicate/O ss1sp $(spname+num+"ovdis")
wave stdis=$(spname+num+"stdis"), st=$(spname+num+"st"), ovdis=$(spname+num+"ovdis")
p=0
do
stdis[p]=sqrt((ss1utmlat[ss1ep[p]]-ss1utmlat[ss1sp[p]])^2+(ss1utmlon[ss1ep[p]]-ss1utmlon[ss1sp[p]])^2)
ovdis[p]=sum(ss1spd, pnt2x(ss1spd,ss1sp[p]), pnt2x(ss1spd,ss1ep[p]))
st[p]=stdis[p]/ovdis[p]
p+=1
while(p<numpnts(ss1sp))
for(p=0; p<numpnts(ss1sp); p+=1)
if(numtype(ss1windspd[p])==2)
stdis[p]=nan; st[p]=nan; ovdis[p]=nan
endif
endfor
//killwaves stdis
end

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////smoothingした経路と比較したときの経路の直線度を計算する関数　smooth_func()//////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
function smooth_func(ss1utmlat, ss1utmlon, ss1sp, ss1ep, ss1spd, spname, num, ss1windspd)   /////   出力wave
string spname, num
wave ss1utmlat, ss1utmlon, ss1sp, ss1ep, ss1spd, ss1windspd
duplicate/O ss1utmlat $(spname+num+"smthlat"); duplicate/O ss1utmlon $(spname+num+"smthlon")
wave smthlat=$(spname+num+"smthlat"), smthlon=$(spname+num+"smthlon")
smooth 100, smthlat, smthlon
	calcspddirturn(smthlat, smthlon)
rename $(spname+num+"smthlon_spd") $(spname+num+"smthspd")
killwaves $(spname+num+"smthlon_turn"), $(spname+num+"smthlon_dir")
duplicate/O ss1sp $(spname+num+"smthst")
wave smthst=$(spname+num+"smthst"), smthspd=$(spname+num+"smthspd")
smthspd = smthspd * 1000/3600
variable p
p=0
do
smthst[p] = sum(smthspd, pnt2x(smthspd, ss1sp[p]), pnt2x(smthspd, ss1ep[p])) / sum(ss1spd, pnt2x(ss1spd, ss1sp[p]), pnt2x(ss1spd, ss1ep[p]))
p+=1
while(p<numpnts(ss1sp))
killwaves smthspd
for(p=0; p<numpnts(ss1windspd); p+=1)
if(numtype(ss1windspd[p])==2)
smthst[p]=nan
endif
endfor
end

//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////      加速度は系から羽ばたきを計算する関数　flap_count()     //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
function flap_count()
wave lay2sp, lay2ep, lay2zw
variable p
duplicate/O lay2zw lay2zw_fil
make/O/D/N=0 coefs; filterfir/DIM=0/HI={0.125, 0.125, 101}/coef coefs, lay2zw_fil
makemask(lay2zw_fil, low=-inf, high=-0.5, mname="ref")
masktopoint(ref, "lay2flapmask")
killwaves ref, lay2zw_fil
make/N=(numpnts(lay2sp))/O lay2flap
p=0
do
	lay2flap[p] = sum(lay2flapmask, pnt2x(lay2flapmask, lay2sp[p]*20), pnt2x(lay2flapmask, lay2ep[p]*20))
p+=1
while(p<numpnts(lay2sp))
end


////////////////////////////////////////////////////////////////////////////////////////////////
//////////          それぞれの区間の鳥の進行方向と風とのなす角度を計算する関数　birdwindangle         //////////
////////////////////////////////////////////////////////////////////////////////////////////////
function birdwindangle()   /////   出力wave head, windang
wave accutmlat,accutmlon,accsp,accep
variable p
duplicate/O accsp acchead; acchead = NAN
p=0
do
acchead[p] = atan2((accutmlat[accep[p]]-accutmlat[accsp[p]]),(accutmlon[accep[p]]-accutmlon[accsp[p]]))
if(acchead[p] > pi/2 && acchead[p] <=pi)
acchead[p] = acchead[p] - 2*pi
endif
p+=1
while (p < numpnts(accsp))
wave accradian, acchead
variable q
duplicate/O accradian accr; accr = accradian + 3/2*pi; duplicate/O acchead acch; acch = acchead + 3/2*pi; duplicate/O accr accwindang
q=0
do
if(accr[q]-acch[q] >= 0 && accr[q]-acch[q] <=pi)
accwindang[q] = accr[q]-acch[q]
endif
if(accr[q]-acch[q] < 0 && accr[q]-acch[q] >=-pi)
accwindang[q] = acch[q]-accr[q]
endif
if(accr[q]-acch[q] > pi && accr[q]-acch[q] <=2*pi)
accwindang[q] = 2*pi-(accr[q]-acch[q])
endif
if(accr[q]-acch[q] < -pi && accr[q]-acch[q] >=-2*pi)
accwindang[q] = 2*pi+(accr[q]-acch[q])
endif
q+=1
while(q<numpnts(accr))
killwaves accr acch
end

/////////////////////////////////////////////////////////////////////////////////////////
//////////          5分区間の飛行経路の蛇行の周期性とその振幅を計算する関数　pathfft         //////////
/////////////////////////////////////////////////////////////////////////////////////////
function pathfft()   /////   出力wave peak, mag
wave ss37spd,ss37sp,ss37ep
make/N=(numpnts(ss37sp)) ss37peak; make/N=(numpnts(ss37sp)) ss37mag
variable p,q
q=0
do
duplicate/R=[ss37sp[q],ss37ep[q]] ss37spd ss37spdcut
duplicate/R=[0,99] ss37spdcut cut
FFT/OUT=3/WINF=Hanning/DEST=FFTref cut
killwaves cut
p=50
do   /////   加算平均
duplicate/R=[p,p+99] ss37spdcut cut
FFT/OUT=3/WINF=Hanning/DEST=cutFFT cut
FFTref = FFTref + cutFFT
killwaves cutFFT cut
p+=50
while(p<=200)
SetScale/P x 0,0.01,"", FFTref
findpeak/Q FFTref
ss37peak[q] = V_peakloc; ss37mag[q] = V_peakval
killwaves ss37spdcut, FFTref
q+=1
while(q<numpnts(ss37sp))
end

///////////////////////////////////////////////////////////////////////////////
//////////          5分区間の羽ばたき時間を加速度から求める関数 flapnum()         //////////
///////////////////////////////////////////////////////////////////////////////
function flapnum()
wave ss37sp, ss37ep, ss37flapmask
variable p
duplicate ss37sp ss37flap
p=0
do
//duplicate/R=[ss37sp[p]*20,ss37ep[p]*20] ss37flapmask flap1	//SSの場合
//ss37flap[p] = sum(flap1)/20	//SSの場合
duplicate/R=[ss37sp[p],ss37ep[p]] ss37flapmask flap1	//Layの場合
ss37flap[p] = sum(flap1)	//Layの場合
killwaves flap1
p+=1
while(p<numpnts(ss37sp))
end


///////////////////////////////////////////////////////////////////////////////
//////////          Wilson et al. 2007 のAIIを計算する関数　AII()         //////////
///////////////////////////////////////////////////////////////////////////////
function aii(ss37utmlat,ss37utmlon)//i, ID, wavelistを個体ごとに指定する
wave ss37utmlat, ss37utmlon
variable scale, p, D, L, q, i, index
string thewavename, thewavename2, ID
//Wilson et al. 2007 Deep Sea Research参照, m=2
i = 120//maximum scaleの指定
ID = "ss37"//個体ごとのwaveを指定

scale = 2
do

thewavename = ID + "k" + num2istr(scale)
make/N=(numpnts(ss37utmlat)) $thewavename

wave referencewave = $thewavename//個体ごとのwave名を使う

p = scale/2//繰り返し変数
q = scale/2//繰り返しなし変数
do
D = sqrt((ss37utmlat[p-q]-ss37utmlat[p+q])^2+(ss37utmlon[p-q]-ss37utmlon[p+q])^2)//scale=2のときの直線移動距離
L = sqrt((ss37utmlat[p-q]-ss37utmlat[p])^2+(ss37utmlon[p-q]-ss37utmlon[p])^2)+sqrt((ss37utmlat[p]-ss37utmlat[p+q])^2+(ss37utmlon[p]-ss37utmlon[p+q])^2)//scale=2のときの水平移動距離
referencewave[p] = D/L//scale=2のときの蛇行度
p+=1
while(p<numpnts(ss37utmlat)-scale/2)

NAN2Zero(referencewave)//Ethographerのコマンド
scale += 2
while(scale<=i)

thewavename2 = ID + "AII"
string list = wavelist("ss37k*",";","")//個体ごとのwaveを指定
concatenate list, $thewavename2
edit $thewavename2

for(index=0; index<i/2; index+=1)//行列を作るのに用いたwaveを消す
killwaves $(stringfromlist(index,list))
endfor

setscale/P x leftx(ss37utmlat),deltax(ss37utmlat),"dat", $thewavename2//x軸のスケールをそろえる

end



///////////////////////////////////////////////////////////////////////////
//////////          Gipsyから緯度経度と時間を読み込む関数　timesub()         //////////
///////////////////////////////////////////////////////////////////////////
function timesub()
wave ss37time, ss37lat, ss37lon
variable i, p
make/N=(numpnts(ss37time)) ss37timesub
for(i=0;i<numpnts(ss37time);i+=1)
ss37timesub[i]=ss37time[i+1]-ss37time[i]
endfor
p=0
do
if(ss37timesub[p]==1||ss37timesub[p]==-86399)
p+=1
elseif(ss37timesub[p]<1)
deletepoints p+1,-ss37timesub[p]+1,ss37time, ss37lat, ss37lon, ss37timesub
p+=1
elseif(ss37timesub[p]>1)
insertpoints p+1,ss37timesub[p]-1,ss37time, ss37lat, ss37lon, ss37timesub
p+=ss37timesub[p]
endif
while(p<numpnts(ss37time)-1)
end

//////////////////////////////////////////////////////////////////////////////////
//////////          GPSデータと風向風速点の時間データを作成する関数　instime()         //////////
//////////////////////////////////////////////////////////////////////////////////
function instime(ss37time, ss37lat, ss37sp, ss37ep)//最初の1点はあらかじめleftxで確認し、*timeに入力しておくこと
wave ss37time, ss37lat, ss37sp, ss37ep
variable p
string thewavename, ID = "ss37"//個体ごとに値を入力
for(p=1;p<numpnts(ss37lat)-1;p+=1)
insertpoints p,1,ss37time
ss37time[p]=ss37time[p-1]+1
endfor
thewavename = ID + "windtime"
make/N=(numpnts(ss37sp)) $thewavename
wave refwave = $thewavename
for(p=0;p<numpnts(ss37sp);p+=1)
duplicate/R = [ss37sp[p],ss37ep[p]] ss37time timecut
refwave[p] = mean(timecut)
killwaves timecut
endfor
end

////////////////////////////////////////////////////////////////////////////////////////////
//////////          　進行方向と風のなす角のwindvectorデータを作成する関数　windheadvec()         //////////
////////////////////////////////////////////////////////////////////////////////////////////
function windheadvec()
wave ss37windang, ss37windspeed //////////　　　進行方向と風のなす角のwindvectorデータを作成する
string thewavename = "ss37" + "windheadvec"//個体ごとに値を入力
duplicate ss37windang a; duplicate ss37windspeed b
a = ss37windang*pi/180 * (-1) +pi/2
b = ss37windspeed * 5
insertpoints numpnts(a), numpnts(a), b
insertpoints 0, numpnts(a), a
make/N = (numpnts(ss37windang),2)/D $thewavename
wave refwave = $thewavename
refwave = a + b
killwaves a b
End


////////////////////////////////////////////////////////////////////////////////////////////
//////////          　速度変化が著しいところを抽出する関数　DSmask()         //////////
////////////////////////////////////////////////////////////////////////////////////////////
function DSmask()//作成されるwave; ssdsmask
wave ssspd
variable p, q
differentiate ssspd/D=ssspd_dif
differentiate ssspd_dif/D=ssspd_dif_dif
duplicate/O ssspd ssspdds; ssspdds = 0
duplicate/O ssspd ssspdds_start; ssspdds_start= nan
duplicate/O ssspd ssspdds_end; ssspdds_end = nan
for(p=1; p<numpnts(ssspd)-1; p+=1)//y=0を超える前の点と超えた後の点を、それぞれ終了点と開始点にするウェーブを作る
	if(ssspd_dif[p-1]<0 && ssspd_dif[p]>=0)
	ssspdds_start[p] = p
	elseif(ssspd_dif[p]<0 && ssspd_dif[p+1]>=0)
	ssspdds_end[p] = p
	else
	endif
endfor
p = 0//ssspdds_start, ssspdds_endのウェーブを調整する（nanを削除し、開始点と終了点をそろえる）
do
	if(numtype(ssspdds_start[p]) == 2)
	deletepoints p,1,ssspdds_start
	else
	p+=1
	endif
while(p<numpnts(ssspdds_start))
p = 0
do
	if(numtype(ssspdds_end[p]) == 2)
	deletepoints p,1,ssspdds_end
	else
	p+=1
	endif
while(p<numpnts(ssspdds_end))
if(ssspdds_end[0]<ssspdds_start[0])
deletepoints 0,1,ssspdds_end
endif
if(numpnts(ssspdds_start)>numpnts(ssspdds_end))
deletepoints numpnts(ssspdds_start)-1,1,ssspdds_start
endif
q = 0//対地速度10m/s以上かつ加速度1m/s以上のピークが見られるサイクルを抽出する
variable a, b
do
	if(wavemax(ssspd, pnt2x(ssspd, ssspdds_start[q]), pnt2x(ssspd, ssspdds_end[q])) >= 10 && wavemax(ssspd_dif, pnt2x(ssspd_dif, ssspdds_start[q]), pnt2x(ssspd_dif, ssspdds_end[q])) >= 1)
	wavestats/Q/R=[ssspdds_start[q], ssspdds_end[q]] ssspd_dif
	a = v_maxrowloc
	b = 0; p = a+1
	do
	if(ssspd_dif_dif[p] > 0)
	b = p
	endif
	p+=1
	while(b<=0)	
	ssspdds[a, b] = 1
	endif
	q+=1
while(q<numpnts(ssspdds_start))
resample/rate=4 ssspdds//マスクデータを4Hzにするかどうかを決める
MakeMask(ssspdds, Low=0.5, High=inf, mName="ss_4hzdsmask")//抽出したサイクルを示すマスクを作る
killwaves ssspdds_start, ssspdds_end, ssspdds
end

////////////////////////////////////////////////////////////////////////////////////////////
//////////          　ダイナミックソアリングしているときとしていないときの羽ばたき数の違い　flaprate()         /////////
////////////////////////////////////////////////////////////////////////////////////////////
function flaprate()
wave accflapmask, accdsmask, accflight, acc_analysismask//acc_analysismaskはcommuting flightの区間だけを抽出したもの
selectdatabymask("accdsmask", "acc_analysismask")
selectdatabymask("accflapmask", "acc_analysismask")
selectdatabymask("accflight", "acc_analysismask")
makemask(accdsmask_select, low=0.5, high=inf, mname="accdsmask1"); killwaves accdsmask_select//commuting flight の区間だけのデータを抽出する
makemask(accflapmask_select, low=0.5, high=inf, mname="accflapmask1"); killwaves accflapmask_select
makemask(accflight_select, low=0.5, high=inf, mname="accflight1"); killwaves accflight_select
duplicate accdsmask1 accdsmask1_inv
maskinverse(accdsmask1_inv)
selectdatabymask("accdsmask1", "accflight1")//accdsmask1_selectができる
selectdatabymask("accflapmask1", "accflight1")//accflapmask1_selectができる
selectdatabymask("accdsmask1_inv", "accflight1")//accdsmask1_inv_selectができる
makemask(accdsmask1_select, low=0.5, high=inf, mname="accdsmask2"); killwaves accdsmask1_select//accdsmask2というflight区間だけのDSのマスクができる
selectdatabymask("accflapmask1_select", "accdsmask2")//accflapmask1_select_selectができる。flight区間のDS区間の羽ばたき
nan2zero(accflapmask1_select_select)
print sum(accflapmask1_select_select), "flight区間のDS区間の羽ばたきの数"
print sum(accdsmask2), "flight区間のDS区間の総数"
print sum(accflapmask1_select_select)/sum(accdsmask2), "flight区間内のDS区間の羽ばたき割合"
killwaves accdsmask2, accflapmask1_select_select
makemask(accdsmask1_inv_select, low=0.5, high=inf, mname="accdsmask1_inv2"); killwaves accdsmask1_inv_select//accdsmask1_inv2というflight区間のDS以外の部分のマスクを作る
selectdatabymask("accflapmask1_select", "accdsmask1_inv2")//accflapmask1_select_selectができる。flight区間のDS以外の部分の羽ばたき
nan2zero(accflapmask1_select_select)
print sum(accflapmask1_select_select), "flight区間のDS以外の区間の羽ばたき数"
print sum(accdsmask1_inv2), "flight区間のDS以外の総数"
print sum(accflapmask1_select_select)/sum(accdsmask1_inv2), "flight区間内のDS以外の区間の羽ばたき割合"
killwaves accdsmask1_inv2, accflapmask1_select_select, accflapmask1_select, accdsmask1_inv, accdsmask1, accflapmask1, accflight1
end


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          移動区間を1分間隔で区切り、DSの有無や羽ばたき割合を計算する関数 flappingcount()         //////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
function flappingcount() //commutesp, commuteep, commuteds, commutedsflap, commutenondsflap, commuteflap
///移動区間の開始点と終了点を求める
wave lay1_4hzcommutemask
variable p, n, s, q
	StatMask(lay1_4hzcommutemask, lay1_4hzcommutemask)
	GetStatData(Stat_Result, "Stat", "startRow")
	GetStatData(Stat_Result, "Stat", "endRow")
	wave startrow_stat, endrow_stat
		make/O/N=3000 lay1comsp, lay1comep; lay1comsp = nan; lay1comep = nan
		p = 0
		n = 0
		do
		q = 0
		do
		lay1comsp[p] = startrow_stat[n] + 60*4*q
		lay1comep[p] = startrow_stat[n] + 60*4*q + 60*4
		p+=1; q+=1
		while(q<((endrow_stat[n]-startrow_stat[n])-60*4)/(60*4))
		n+=1
		while(n<numpnts(startrow_stat))
			s=0
			do
			if(numtype(lay1comsp[S])==2)
			deletepoints S, 1, lay1comsp, lay1comep
			else
			s+=1
			endif
			while(s<numpnts(lay1comsp))
	killwaves stat_result, startrow_stat, endrow_stat
///各区間の羽ばたきやDSの時間割合を求める
wave lay1_4hzdsmask, lay1_4hzflapmask
variable a, b, c, d
duplicate/O lay1comsp, lay1comds, lay1comdsflap, lay1comnondsflap, lay1comflap
	for(p=0; p<numpnts(lay1comsp); p+=1)
	duplicate/O/R=[lay1comsp[p], lay1comep[p]] lay1_4hzdsmask dsmask1
	duplicate/O/R=[lay1comsp[p], lay1comep[p]] lay1_4hzflapmask flapmask1
		a=0; b=0; c=0; d=0
		for(n=0; n<numpnts(dsmask1); n+=1)
		if(dsmask1[n]==1 && flapmask1[n]==1)
		a = a+1
		elseif(dsmask1[n]==1 && flapmask1[n]==0)
		b = b+1
		elseif(dsmask1[n]==0 && flapmask1[n]==1)
		c = c+1
		elseif(dsmask1[n]==0 && flapmask1[n]==0)
		d = d+1
		endif
		endfor
	lay1comds[p] = sum(dsmask1)/numpnts(dsmask1)//区間内のDSしている時間割合
	lay1comdsflap[p] = a/(a+b)//DSしているときの羽ばたき時間割合
	lay1comnondsflap[p] = c/(c+d)//DSしていないときの羽ばたき時間割合
	lay1comflap[p] = (a+c)/(a+b+c+d)//区間内の羽ばたき時間割合
	endfor
	killwaves dsmask1, flapmask1
nan2zero(lay1comdsflap)
edit lay1comsp lay1comep lay1comds lay1comdsflap lay1comnondsflap lay1comflap
end

////////////////////////////////////////////////////////////////////////////////////////////////
//////////          　ダイナミックソアリングをしているときの羽ばたき量の比較をする関数　dsflap()         /////////
////////////////////////////////////////////////////////////////////////////////////////////////
function dsflap() //新しいwaveはできない
wave acccommutemask, accflapmask, accdsmask
variable p,a,b,c,d
a=0; b=0; c=0; d=0
for(p=0; p<numpnts(accflapmask); p+=1)
if(acccommutemask[p]==1 && accdsmask[p]==1 && accflapmask[p]==1)//ダイナミックソアリングしていて羽ばたくとき
a=a+1
elseif(acccommutemask[p]==1 && accdsmask[p]==1 && accflapmask[p]==0)//ダイナミックソアリングしていて羽ばたかないとき
b=b+1
elseif(acccommutemask[p]==1 && accdsmask[p]==0 && accflapmask[p]==1)//ダイナミックソアリングしていなくて羽ばたくとき
c=c+1
elseif(acccommutemask[p]==1 && accdsmask[p]==0 && accflapmask[p]==0)//ダイナミックソアリングしていなくて羽ばたかないとき
d=d+1
endif
endfor
print a, b, c, d
print "DS中の羽ばたき割合", a/(a+b)
print "DS以外の羽ばたき割合", c/(c+d)
print "移動飛行中の羽ばたき割合", (a+c)/(a+b+c+d)
end

////////////////////////////////////////////////////////////////////////////////////////////////
//////////          　ダイナミックソアリングをしているときの羽ばたき量の比較をする関数　dsflap()         /////////
////////////////////////////////////////////////////////////////////////////////////////////////
function windresamp()
wave accwindspeed, acctime, acclat, accradian
variable p
duplicate acclat accwindspeed_resamp; accwindspeed_resamp = 0
duplicate acclat accradian_resamp; accradian_resamp = 0
for(p=0; p<numpnts(accwindspeed); p+=1)
accwindspeed_resamp[x2pnt(acclat, acctime[p]-14*4), x2pnt(acclat, acctime[p]+15*4)] = accwindspeed[p]
accradian_resamp[x2pnt(acclat, acctime[p]-14*4), x2pnt(acclat, acctime[p]+15*4)] = accradian[p]
endfor
end

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

Window SS_mag_vs_windang() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(43,73.5,438,281.5) ss1mag vs ss1windang as "SS mag vs windang all"
	AppendToGraph ss2mag vs ss2windang
	AppendToGraph ss3mag vs ss3windang
	AppendToGraph ss4mag vs ss4windang
	AppendToGraph ss5mag vs ss5windang
	AppendToGraph ss6mag vs ss6windang
	AppendToGraph ss7mag vs ss7windang
	AppendToGraph ss11mag vs ss11windang
	AppendToGraph ss12mag vs ss12windang
	AppendToGraph ss13mag vs ss13windang
	AppendToGraph ss14mag vs ss14windang
	AppendToGraph ss15mag vs ss15windang
	AppendToGraph ss16mag vs ss16windang
	ModifyGraph mode=3
	ModifyGraph marker=19
	ModifyGraph rgb(ss1mag)=(0,0,0)
	ModifyGraph msize=2
	ModifyGraph zColor(ss1mag)={ss1windspeed,0,10,Rainbow,1},zColor(ss2mag)={ss2windspeed,0,10,rainbow,1}
	ModifyGraph zColor(ss3mag)={ss3windspeed,0,10,rainbow,1},zColor(ss4mag)={ss4windspeed,0,10,rainbow,1}
	ModifyGraph zColor(ss5mag)={ss5windspeed,0,10,rainbow,1},zColor(ss6mag)={ss6windspeed,0,10,rainbow,1}
	ModifyGraph zColor(ss7mag)={ss7windspeed,0,10,rainbow,1},zColor(ss11mag)={ss11windspeed,0,10,rainbow,1}
	ModifyGraph zColor(ss12mag)={ss12windspeed,0,10,rainbow,1},zColor(ss13mag)={ss13windspeed,0,10,rainbow,1}
	ModifyGraph zColor(ss14mag)={ss14windspeed,0,10,rainbow,1},zColor(ss15mag)={ss15windspeed,0,10,rainbow,1}
	ModifyGraph zColor(ss16mag)={ss16windspeed,0,10,rainbow,1}
	ModifyGraph mirror=2
	ModifyGraph font="Arial"
	ModifyGraph fSize=10
	ModifyGraph standoff=0
	ModifyGraph manTick(bottom)={0,90,0,0},manMinor(bottom)={0,0}
	Label left "\\F'arial'\\Z10Magnitude"
	Label bottom "\\F'arial'\\Z10Wind direction (degrees)"
	SetAxis left 0,*
	SetAxis bottom 0,360
	ColorScale/C/N=text0/F=0/B=1/A=MC/X=34.80/Y=34.19 trace=ss1mag, vert=0
	ColorScale/C/N=text0 heightPct=5, width=50, tickLen=3, nticks=2, font="Arial"
	ColorScale/C/N=text0 fsize=10, lblMargin=0
	AppendText "\\F'Arial'\\Z10Wind speed (m/s)"
EndMacro

Window Graph0() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(35,42.5,746,452) ssaccspd
	AppendToGraph/R=ds ssaccspdds
	AppendToGraph ssaccspd_fil
	ModifyGraph mode(ssaccspd)=4
	ModifyGraph marker(ssaccspd)=19
	ModifyGraph zmrkSize(ssaccspd)={ssaccflapmask,*,*,0,1.5}
	ModifyGraph zColor(ssaccspd)={ssaccflapmask,*,*,RedWhiteBlue,1}
	ModifyGraph grid(ds)=1
	ModifyGraph zero(left)=2,zero(ds)=1
	ModifyGraph font(left)="Arial",font(ds)="Arial"
	ModifyGraph sep(ds)=2
	ModifyGraph fSize(left)=15,fSize(ds)=15
	ModifyGraph standoff=0
	ModifyGraph gridHair(ds)=1
	ModifyGraph lblPos(ds)=51
	ModifyGraph lblLatPos(ds)=8
	ModifyGraph freePos(ds)=0
	ModifyGraph axisEnab(left)={0,0.5}
	ModifyGraph axisEnab(ds)={0.5,1}
	ModifyGraph manTick(ds)={0,10,0,0},manMinor(ds)={1,0}
	ModifyGraph dateInfo(bottom)={0,0,0}
	Label left "\\F'Arial'\\Z15Ground velocity (m/s)"
	Label ds "\\F'Arial'\\Z15Sum of absolute speed"
	SetAxis bottom 3460343772.14094,3460344056.90571
	Cursor/P A ssaccspd 22955
	ShowInfo
EndMacro


////////////////////////////////////////////////////////////////////////////////////////////////
//////////          それぞれの区間の鳥の進行方向と風とのなす角度を計算する関数　birdwindangle         //////////
////////////////////////////////////////////////////////////////////////////////////////////////
function birdwindangle2()   /////   出力wave head, windang windang180
wave comsp, comep, ssaccutmlat, ssaccutmlon
variable p
duplicate/O comsp comhead; comhead = NAN
	p=0
	do
	comhead[p] = atan2((ssaccutmlat[comep[p]/4]-ssaccutmlat[comsp[p]/4]),(ssaccutmlon[comep[p]/4]-ssaccutmlon[comsp[p]/4]))
	if(comhead[p] > pi/2 && comhead[p] <=pi)
	comhead[p] = comhead[p] - 2*pi
	endif
	p+=1
	while (p < numpnts(comsp))
wave comradian, comhead
variable q
duplicate/O comradian comr; comr = comradian + 3/2*pi; duplicate/O comhead comh; comh = comhead + 3/2*pi; duplicate/O comr comwindang
	q=0
	do
	if(comr[q]-comh[q] >= 0 && comr[q]-comh[q] <=pi)
	comwindang[q] = comr[q]-comh[q]
	endif
	if(comr[q]-comh[q] < 0 && comr[q]-comh[q] >=-pi)
	comwindang[q] = comh[q]-comr[q]
	endif
	if(comr[q]-comh[q] > pi && comr[q]-comh[q] <=2*pi)
	comwindang[q] = 2*pi-(comr[q]-comh[q])
	endif
	if(comr[q]-comh[q] < -pi && comr[q]-comh[q] >=-2*pi)
	comwindang[q] = 2*pi+(comr[q]-comh[q])
	endif
	q+=1
	while(q<numpnts(comr))
killwaves comr comh
duplicate/O comwindang comwindang180
for(p=0; p<numpnts(comwindang); p+=1)
if(comwindang[p] > pi)
comwindang180[p] = 2*pi - comwindang[p]
endif
endfor
end

/////////////////////////////////////////////////////////////////////////
//////////          それぞれの区間の移動距離を計算する関数　move()         //////////
/////////////////////////////////////////////////////////////////////////

function move() // 出力ウェーブ sscomdis
wave sscomsp, sscomep, ssutmlat, ssutmlon, ssspd
variable p, sp1, ep1
make/O/N=(numpnts(sscomsp)) sscomdis, sscomstdis
for(p=0; p < numpnts(sscomsp); p+=1)
sp1 = sscomsp[p]/4
ep1 = sscomep[p]/4
sscomstdis[p] = sqrt((ssutmlat[sp1] - ssutmlat[ep1])^2 + (ssutmlon[sp1] - ssutmlon[ep1])^2)
sscomdis[p] = sum(ssspd, pnt2x(ssspd, sp1), pnt2x(ssspd, ep1))
endfor
end


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////          移動区間を1分間隔で区切り、DSの有無や羽ばたき割合を計算する関数 flappingcount()         //////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
function flappingcount2() //commutesp, commuteep, commuteds, commutedsflap, commutenondsflap, commuteflap
///移動区間の開始点と終了点を求める
wave lay1_4hzcommutemask
variable p, n, s, q
	StatMask(lay1_4hzcommutemask, lay1_4hzcommutemask)
	GetStatData(Stat_Result, "Stat", "startRow")
	GetStatData(Stat_Result, "Stat", "endRow")
	wave startrow_stat, endrow_stat
		make/O/N=3000 lay1comsp, lay1comep; lay1comsp = nan; lay1comep = nan
		p = 0
		n = 0
		do
		q = 0
		do
		lay1comsp[p] = startrow_stat[n] + 60*4*q
		lay1comep[p] = startrow_stat[n] + 60*4*q + 60*4
		p+=1; q+=1
		while(q<((endrow_stat[n]-startrow_stat[n])-60*4)/(60*4))
		n+=1
		while(n<numpnts(startrow_stat))
			s=0
			do
			if(numtype(lay1comsp[S])==2)
			deletepoints S, 1, lay1comsp, lay1comep
			else
			s+=1
			endif
			while(s<numpnts(lay1comsp))
	killwaves stat_result, startrow_stat, endrow_stat
///各区間の羽ばたきやDSの時間割合を求める
wave lay1_4hzdsmask, lay1_4hzflapmask
variable a, b, c, d
duplicate/O lay1comsp, lay1comds, lay1comdsflap, lay1comnondsflap, lay1comflap
	for(p=0; p<numpnts(lay1comsp); p+=1)
	duplicate/O/R=[lay1comsp[p], lay1comep[p]] lay1_4hzdsmask dsmask1
	duplicate/O/R=[lay1comsp[p], lay1comep[p]] lay1_4hzflapmask flapmask1
		a=0; b=0; c=0; d=0
		for(n=0; n<numpnts(dsmask1); n+=1)
		if(dsmask1[n]==1 && flapmask1[n]==1)
		a = a+1
		elseif(dsmask1[n]==1 && flapmask1[n]==0)
		b = b+1
		elseif(dsmask1[n]==0 && flapmask1[n]==1)
		c = c+1
		elseif(dsmask1[n]==0 && flapmask1[n]==0)
		d = d+1
		endif
		endfor
	lay1comds[p] = sum(dsmask1)/numpnts(dsmask1)//区間内のDSしている時間割合
	lay1comdsflap[p] = a/(a+b)//DSしているときの羽ばたき時間割合
	lay1comnondsflap[p] = c/(c+d)//DSしていないときの羽ばたき時間割合
	lay1comflap[p] = (a+c)/(a+b+c+d)//区間内の羽ばたき時間割合
	endfor
	killwaves dsmask1, flapmask1
nan2zero(lay1comdsflap)
edit lay1comsp lay1comep lay1comds lay1comdsflap lay1comnondsflap lay1comflap
end

////////////////////////////////////////////////////////////////
//////////          横風成分と追い風成分を作る　swctwc()         //////////
////////////////////////////////////////////////////////////////
function swctwc(ss31windspd, ss31bwang, spname, num)
wave ss31windspd, ss31bwang
string spname, num
make/O/N=(numpnts(ss31windspd)) $(spname+num+"swc"), $(spname+num+"twc")
wave swc=$(spname+num+"swc"), twc=$(spname+num+"twc")
swc = ss31windspd * sin(ss31bwang*pi/180)
twc = ss31windspd * cos(ss31bwang*pi/180)
end


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

function flap()
wave lay1zwflap, lay1sp, lay1ep
variable p
duplicate/O lay1sp lay1flap
for(p=0; p<numpnts(lay1sp); p+=1)
lay1flap[p]=sum(lay1zwflap, pnt2x(lay1zwflap, lay1sp[p]*20), pnt2x(lay1zwflap, lay1ep[p]*20))/(300*20)
endfor
end

function delete()
wave lay1flap, lay1windspd, lay1bwang, lay1st, lay1swc, lay1twc, lay1id, lay1flapnum
variable p
p=0
do
if(numtype(lay1windspd)==2)
deletepoints p,1,lay1flap, lay1windspd, lay1bwang, lay1st, lay1swc, lay1twc, lay1id, lay1flapnum
else
p+=1
endif
while(p<numpnts(lay1windspd))
end

function temporally()
variable p
p=3
do
wave spd=$("wanspd"+num2str(p)), rad=$("wanrad"+num2str(p))
display spd vs rad
modifygraph mode=3, marker=19, msize=1.5, rgb=(0,0,0), width=8*30, height={aspect,1},font='arial',fsize=12
label bottom "\\F'arial'\\Z12Flight direction (deg)"
label left "\\F'arial'\\Z12Ground speed (m/s)"
setaxis bottom 0,pi*2
setaxis left 0,40
K0 = 15.4;K2 = 1;CurveFit/X=1/H="1010"/NTHR=0 sin  spd /X=rad /D /F={0.950000, 5}
string fit="fit_wanspd"+num2str(p), uc="uc_wanspd"+num2str(p), lc="lc_wanspd"+num2str(p), spdwave="wanspd"+num2str(p)
ModifyGraph lsize($(fit))=2,mode($(UC))=7,hbFill($(UC))=2
ModifyGraph rgb($(UC))=(47872,47872,47872),mode($(LC))=7
ModifyGraph hbFill($(LC))=2,rgb($(LC))=(65535,65535,65535)
ReorderTraces $(spdwave),{$(UC),$(LC)}
string name="wanspd"+num2str(p)+"_wanrad"+num2str(p)
DoWindow/C/T name, name
wave fitwave=$("fit_wanspd"+num2str(p)), wanwindspd, wanwindrad
wavestats fitwave; wanwindspd[p-1] = (v_max-v_min) /2; wanwindrad[p-1] = pi/2 - v_maxloc
p+=1
while(p<=9)
end

function temp_windmake()
wave start2, end2
variable p,q,r
make/O/N=3000 sp2, ep2
p=0;q=0
do
r=0
do
sp2[p]=start2[q]+60+r*300
ep2[p]=start2[q]+60+r*300+300
p+=1
r+=1
while(r<((end2[q]-start2[q])-120-300)/300)
q+=1
while(q<numpnts(start2))
end

function tapioka()
wave sp2, ep2
variable p
make/N=(numpnts(sp2))/O/D windrad2=nan, windspd2=nan, lonave2=nan, latave2=nan, sinaic2=nan, linaic2=nan, sinlin2=nan
p=0
do
	duplicate/O/R = [sp2[p],ep2[p]] lat2 latcut
	duplicate/O/R = [sp2[p],ep2[p]] lon2 loncut
	lonave2[p] = mean(loncut); latave2[p] = mean(latcut)
	duplicate/O/R = [sp2[p],ep2[p]] spd2 speed
	duplicate/O/R = [sp2[p],ep2[p]] rad2 angle
		display speed vs angle; setaxis bottom 0, 2*pi
		K0=10;K2 = 1;//////////     choose air speed
		CurveFit/Q=1/L=360 /X=1 /H="1010"/NTHR=0 sin speed /X=angle /D /R //////////     sine fitting
		WaveStats/Q fit_speed
		windrad2[p] = pi/2 - V_maxloc
		windspd2[p] = (V_max - V_min)/2
		wave Res_speed
		duplicate/O Res_speed res; res = Res_speed ^ 2
		sinaic2[p] = numpnts(res)*((log(2*pi*sum(res)/numpnts(res))/log(e))+1)+2*(1+2)
		dowindow/K graph1 //////////     choose graph number
	killwaves W_coef W_sigma fit_speed res_speed res
		display speed vs angle; setaxis bottom 0, 2*pi
		K1 = 0;
		CurveFit/Q=1/L=360 /X=1/H="01"/NTHR=0 line  speed /X=angle /D /R //////////     line fitting
		wave Res_speed
		duplicate/O Res_speed res
		res = Res_speed ^ 2
		linaic2[p] = numpnts(res)*((log(2*pi*sum(res)/numpnts(res))/log(e))+1)+2*(0+2)
	if(sinaic2[p]>=linaic2[p]-2) //////////     compare sine and line
	sinlin2[p] = 1
	else
	sinlin2[p] = 0
	endif
		dowindow/K graph1 //////////     choose graph number
	if(windspd2[p] > 100 || windspd2[p] < 0.01)
	windrad2[p]=0; windspd2[p]=0; lonave2[p]=0; latave2[p]=0
	sinaic2[p]=0; linaic2[p]=0; sinlin2[p]=1
	endif
killwaves speed angle latcut loncut W_coef fit_speed W_sigma res_speed res
p+=1
while (p < numpnts(sp2))
edit windrad2 windspd2 latave2 lonave2 sp2 ep2
variable s
	s=0
	do
	if(numtype(windrad2[S])==2)
	deletepoints S, 1, windrad2, windspd2, latave2, lonave2, sinaic2, linaic2, sinlin2
	else
	S+=1
	endif
	while(S<numpnts(windrad2))
	S=0
	do
	if(sinlin2[S]==1)
	windrad2[s]=nan; windspd2[s]=nan; latave2[s]=nan; lonave2[s]=nan; sinaic2[s]=nan; linaic2[s]=nan; sinlin2[s]=nan;
	endif
	s+=1
	while(S<numpnts(sinlin2))
	killwaves sinlin2, sinaic2, linaic2
//////////　　　windvectorデータを作成する
duplicate/O windspd2 a; a = windspd2 * 5
concatenate {a, windrad2}, windvec2
killwaves a
End


function pathwindangle2()   /////   出力wave vechead, vecwindang
wave sp2, ep2
variable p
make/O/N=(numpnts(sp2)) head2, length2
p=0
do
	duplicate/O/R=[sp2[p], ep2[p]] rad2 radcut
	duplicate/O/R=[sp2[p], ep2[p]] spd2 spdcut
	concatenate/O {radcut, spdcut}, veccut
	statscircularmeans/Q veccut
	wave w_circularmeans
	head2[p] = w_circularmeans[2]
	length2[p] = w_circularmeans[1]
	if(head2[p]>=0 && head2[p] <pi)//ベクトルの向きを風ベクトル（x軸が0度）と合わせる
	head2[p] = pi/2 - head2[p]
	else
	head2[p] = (5*pi)/2 - head2[p]
	endif 
	p+=1
while(p<numpnts(sp2))
killwaves length2
variable q
wave windrad2
duplicate/O windrad2 r2; r2 = windrad2 + 3/2*pi; duplicate/O head2 h2; h2 = head2 + 3/2*pi
duplicate/O r2 windang2
q=0
do
	if(r2[q]-h2[q] >= 0 && r2[q]-h2[q] <=pi)
	windang2[q] = r2[q]-h2[q]
	endif
	if(r2[q]-h2[q] < 0 && r2[q]-h2[q] >=-pi)
	windang2[q] = h2[q]-r2[q]
	endif
	if(r2[q]-h2[q] > pi && r2[q]-h2[q] <=2*pi)
	windang2[q] = 2*pi-(r2[q]-h2[q])
	endif
	if(r2[q]-h2[q] < -pi && r2[q]-h2[q] >=-2*pi)
	windang2[q] = 2*pi+(r2[q]-h2[q])
	endif
	q+=1
	while(q<numpnts(r2))
killwaves r2 h2
duplicate/O windang2 bwang2
for(p=0; p<numpnts(windang2); p+=1)
if(windang2[p] > pi)
bwang2[p] = 2*pi - windang2[p]
endif
endfor
killwaves windang2
killwaves radcut spdcut veccut w_circularmeans
bwang2 = bwang2 *180/pi
end


function insert()
wave wan4delta, wan4lat,wan4lon,wan4alt
variable p
p=0
do
if(wan4delta[p]>1)
insertpoints p,(wan4delta[p]-1),wan4lat,wan4lon,wan4alt
p=p+wan4delta[p]
else
p=p+1
endif
while(p<sum(wan4delta))
zero2nan(wan4lat)
zero2nan(wan4lon)
calcinterp(wan4lat,1)
calcinterp(wan4lon,1)
rename wan4lat wan4latref
rename wan4lon wan4lonref
rename wan4lat_intp wan4lat
rename wan4lon_intp wan4lon
end


function distance()
string id
variable n
n=7
do
id=num2istr(n)
wave lat=$("wan"+id+"utmlat"), lon=$("wan"+id+"utmlon"),spd=$("wan"+id+"windspd"),sp=$("wan"+id+"sp"),ep=$("wan"+id+"ep")
duplicate/O spd $("wan"+id+"dis")
wave dis=$("wan"+id+"dis")
dis=nan
variable lata,latb,lona,lonb,p
p=0
do
lata=lat[sp[p]]
latb=lat[ep[p]]
lona=lon[sp[p]]
lonb=lon[ep[p]]
dis[p]=sqrt((latb-lata)^2+(lonb-lona)^2)
p+=1
while(p<numpnts(spd))
n+=1
while(n<=7)
end



duplicate ss16lat ss16gpstime
ss16gpstime = nan
ss16gpstime = leftx(ss16lat)+p