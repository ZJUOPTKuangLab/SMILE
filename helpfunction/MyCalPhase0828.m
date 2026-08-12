function [phase, amp, offset] = MyCalPhase0828(intlist, phase0)
    int1 = intlist(1);
    int2 = intlist(2);
    int3 = intlist(3);
    cp = cos(phase0);
    sp = sin(phase0);
    
    offset = (int2*cp - (int1+int3)/2)/(cp-1);
    int1 = int1 - offset;
    int2 = int2 - offset;
    int3 = int3 - offset;
    
    amp = sqrt(int2^2 + (int1-int2*cp)^2/(sp^2));
    
    sx = int2/amp;
    cx = -(int1 - int2*cp)/sp/amp;
    
    tp = acos(cx);   % 注意这里acos的范围
%     对于 X 在区间 [-1, 1] 内的实数值，acos(X) 返回区间 [0, π] 内的值。
%     对于 X 在区间 [-1,1] 之外的实数值以及 X 的复数值，acos(X) 返回复数值。
% 所以下面这一步是实现-pi到pi的相位的取值
    if sx <0
        phase = -real(tp);
    else
        phase = real(tp);
    end
end