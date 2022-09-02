function f = GMM_objective(theta,Z,W,residual_fn)

u = residual_fn(theta);
Zu = Z'*u;
f = Zu'*(W\Zu); % same as Zu'*inv(W)*Zu
