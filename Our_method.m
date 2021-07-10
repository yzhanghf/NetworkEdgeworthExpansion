function [EdgeworthCDFValues, EdgeworthCI] = Our_method(A, MotifName, opts)
	
	StdNormal  = opts.StdNormal;
	TestPoints = opts.TestPoints;
	r = opts.r;
	if(isfield(opts, 'alpha')); alpha = opts.alpha; else alpha=-100; end;
	
	n = size(A,1);
	
	% Estimation
	Uhat = Motif(A,0,MotifName);
	% Compute auxiliary quantities;
	G2X1X2 = Motif(A,2,MotifName);
	G1X1 = sum(G2X1X2,2);
	samplesd = sqrt(moment(G1X1,2));
	mu3 = moment(G1X1,3);
	G1X = G1X1*ones(1,n);
	%%% compute E[g1(X1)g1(X2)g2(X1,X2)]
	g1g1g2 = mean(mean(G2X1X2.*G1X.*G1X'));
	
	EdgeworthCDFValues = cdf(StdNormal, TestPoints) + pdf(StdNormal, TestPoints)/sqrt(n)  .*  ...
					(  (1/3*TestPoints.^2+1/6) * (mu3/samplesd^3) + (r-1)/2*(TestPoints.^2+1)* (g1g1g2/samplesd^3)  );
				 
	if(alpha<-1)
		EdgeworthCI = [];
	else
		z_alpha = norminv(alpha/2);
		z_oneminusalpha = norminv(1-alpha/2);
		EdgeworthCI(1) = z_alpha - 1/sqrt(n)  .*  ...
					((1/3*z_alpha^2+1/6) * (mu3/samplesd^3) + (r-1)/2*(z_alpha^2+1)* (g1g1g2/samplesd^3));
		EdgeworthCI(2) = z_oneminusalpha - 1/sqrt(n)  .*  ...
					((1/3*z_oneminusalpha^2+1/6) * (mu3/samplesd^3) + (r-1)/2*(z_oneminusalpha^2+1)* (g1g1g2/samplesd^3));
	end
	
end