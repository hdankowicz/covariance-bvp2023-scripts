function [prob, atlas, cseg, correct] = predict(atlas, prob, cseg) %#ok<INUSD>
%PREDICT   Compute predictor.
%
% Construct curve segment projection condition and theta-method predictor
% and correct.
%
% Identical to atlas2d_v2.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: predict.m 3134 2019-07-13 15:13:42Z hdankowicz $

chart  = atlas.charts{atlas.boundary(1)};
s      = chart.G*chart.s(:,chart.bv(1)); % First available direction
nrms   = norm(s);
h      = atlas.cont.Rmarg*chart.R*nrms;
s      = s/nrms;
prcond = struct('x', chart.x, 'TS', chart.TS, 's', s, 'h', h);
th     = atlas.cont.theta;
if th>=0.5 && th<=1 % A two-step algorithm
  xp           = chart.x + (th*h)*(chart.TS*s);
  [prob, cseg] = CurveSegment.create(prob, chart, prcond, xp);
  [prob, ch2]  = cseg.update_TS(prob, cseg.curr_chart);
  s            = h*ch2.TS'*chart.TS*s;
  h            = norm(s);
  s            = s/h;
  xp           = chart.x + h*ch2.TS*s; % Compute new predictor
  prcond       = struct('x', chart.x, 'TS', ch2.TS, 's', s, 'h', h); % Compute new projection condition
else
  xp           = chart.x + h*chart.TS*s;
end
[prob, cseg] = CurveSegment.create(prob, chart, prcond, xp);
[prob, chart]  = cseg.update_t(prob, cseg.ptlist{1}); % Update tangent vector for event detection
[prob, chart]  = cseg.update_p(prob, chart);          % Update continuation parameters that depend on tangent vector
cseg.ptlist{1} = chart;

correct      = true;

end
