import("music.lib");
import("filter.lib");
import("../modes/roundBigCenterFreq.lib");
import("../modes/roundBigCenterGain.lib");
import("../modes/roundBigCenterT60.lib");
import("../modes/roundMidCenterFreq.lib");
import("../modes/roundMidCenterGain.lib");
import("../modes/roundMidCenterT60.lib");
import("../modes/roundSmallCenterFreq.lib");
import("../modes/roundSmallCenterGain.lib");
import("../modes/roundSmallCenterT60.lib");
import("../modes/semBigCenterFreq.lib");
import("../modes/semBigCenterGain.lib");
import("../modes/semBigCenterT60.lib");
import("../modes/semMidCenterFreq.lib");
import("../modes/semMidCenterGain.lib");
import("../modes/semMidCenterT60.lib");
import("../modes/semSmallCenterFreq.lib");
import("../modes/semSmallCenterGain.lib");
import("../modes/semSmallCenterT60.lib");
import("../modes/squareBigCenterFreq.lib");
import("../modes/squareBigCenterGain.lib");
import("../modes/squareBigCenterT60.lib");
import("../modes/squareMidCenterFreq.lib");
import("../modes/squareMidCenterGain.lib");
import("../modes/squareMidCenterT60.lib");
import("../modes/squareSmallCenterFreq.lib");
import("../modes/squareSmallCenterGain.lib");
import("../modes/squareSmallCenterT60.lib");

N = 36;

x = hslider("x",0,0,1,0.001)*2;
y = hslider("y",0,0,1,0.001)*2;
gate = button("gate");

modeFilter(f,t60) = tf2(b0,b1,b2,a1,a2)
with{
	b0 = 1;
	b1 = 0;
	b2 = -1;
	w = 2*PI*f/SR;
	r = pow(0.001,1/float(t60*SR));
	a1 = -2*r*cos(w);
	a2 = r^2;
};

modeFreq(i) = zxd + (zxu-zxd)*select2(y<1,(y-1),y)
with{
	zx0 = squareSmallCenterFreq(i) + (semSmallCenterFreq(i)-squareSmallCenterFreq(i))*x;
	zx1 = semSmallCenterFreq(i) + (roundSmallCenterFreq(i)-semSmallCenterFreq(i))*(x-1);
	zx2 = squareMidCenterFreq(i) + (semMidCenterFreq(i)-squareMidCenterFreq(i))*x;
	zx3 = semMidCenterFreq(i) + (roundMidCenterFreq(i)-semMidCenterFreq(i))*(x-1);
	zx4 = squareBigCenterFreq(i) + (semBigCenterFreq(i)-squareBigCenterFreq(i))*x;
	zx5 = semBigCenterFreq(i) + (roundBigCenterFreq(i)-semBigCenterFreq(i))*(x-1);
	zxd = select2(y<1,select2(x<1,zx3,zx2),select2(x<1,zx1,zx0));
	zxu = select2(y<1,select2(x<1,zx5,zx4),select2(x<1,zx3,zx2));
};

modeGain(i) = zxd + (zxu-zxd)*select2(y<1,(y-1),y) 
with{
	zx0 = squareSmallCenterGain(i) + (semSmallCenterGain(i)-squareSmallCenterGain(i))*x;
	zx1 = semSmallCenterGain(i) + (roundSmallCenterGain(i)-semSmallCenterGain(i))*(x-1);
	zx2 = squareMidCenterGain(i) + (semMidCenterGain(i)-squareMidCenterGain(i))*x;
	zx3 = semMidCenterGain(i) + (roundMidCenterGain(i)-semMidCenterGain(i))*(x-1);
	zx4 = squareBigCenterGain(i) + (semBigCenterGain(i)-squareBigCenterGain(i))*x;
	zx5 = semBigCenterGain(i) + (roundBigCenterGain(i)-semBigCenterGain(i))*(x-1);
	zxd = select2(y<1,select2(x<1,zx3,zx2),select2(x<1,zx1,zx0));
	zxu = select2(y<1,select2(x<1,zx5,zx4),select2(x<1,zx3,zx2));
};

modeT60(i) = zxd + (zxu-zxd)*select2(y<1,(y-1),y) //: min(0.1)
with{
	zx0 = squareSmallCenterT60(i) + (semSmallCenterT60(i)-squareSmallCenterT60(i))*x;
	zx1 = semSmallCenterT60(i) + (roundSmallCenterT60(i)-semSmallCenterT60(i))*(x-1);
	zx2 = squareMidCenterT60(i) + (semMidCenterT60(i)-squareMidCenterT60(i))*x;
	zx3 = semMidCenterT60(i) + (roundMidCenterT60(i)-semMidCenterT60(i))*(x-1);
	zx4 = squareBigCenterT60(i) + (semBigCenterT60(i)-squareBigCenterT60(i))*x;
	zx5 = semBigCenterT60(i) + (roundBigCenterT60(i)-semBigCenterT60(i))*(x-1);
	zxd = select2(y<1,select2(x<1,zx3,zx2),select2(x<1,zx1,zx0));
	zxu = select2(y<1,select2(x<1,zx5,zx4),select2(x<1,zx3,zx2));
};

mode(i) = modeFilter(modeFreq(i),modeT60(i))*modeGain(i);

body = _ <: par(i,N,mode(i)) :> _;

impulsify = _ <: _,mem : - : >(0);

process = gate : impulsify : body;








