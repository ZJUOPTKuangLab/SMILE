function [lp, resnormx, resultp] = FitPhasePlane3D0828_v2(x,y, z, p2, lp0, fast)
%     parfor m=1:length(aa(:))
%         [lp, resnormx] = FitPhasePlane(xs,ys,p2s, [aa(m) bb(m) 0], 1);
%         resnormmap(m) = resnormx;
%         lplist(m,:) = lp;
%     end

% 其中xs为3552*1的x坐标数据，ys为3552*1的y坐标数据，p2s为3552*1的相应的phase数据
% lp0(1) 为aa(m), lp0(2)为bb(m),m从1到16*22，也就是说把aa还有bb给取遍了
    if nargin<6
        fast = 0;
    end
    
    if fast==0
        options = optimset('Display','off','MaxFunEvals',2000,'MaxIter',1000,'TolFun',1e-10,'LargeScale','on');
    else
        options = optimset('Display','off','MaxFunEvals',500,'MaxIter',100,'TolFun',1e-6,'LargeScale','off'); %fast==1时候
    end
	
    [lpx,resnormx,~,~]=lsqcurvefit(@(xp, xdata)CalPhaseSin3D0828(xp, x, y, z), ...
    [lp0(1) lp0(2) lp0(3) lp0(4)],x,[sin(p2) cos(p2)],[],[],options);

    lp = [lpx(1), lpx(2), lpx(3), checkPhase(lpx(4))];
    if nargout >=3
        resultp = CalPhase(lp, x, y);
    end
end