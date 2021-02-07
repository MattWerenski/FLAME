fbp = [0.409	0.351	0.355	0.371	0.343
0.264	0.222	0.231	0.232	0.223
0.24	0.192	0.194	0.199	0.205
0.402	0.36	0.374	0.383	0.385
0.258	0.227	0.242	0.238	0.247
0.25	0.199	0.221	0.215	0.222
0.48	0.468	0.454	0.474	0.467
0.32	0.319	0.315	0.329	0.317
0.37	0.356	0.336	0.362	0.351];

mbp = [0.396	0.345	0.355	0.361	0.34
0.246	0.209	0.225	0.233	0.209
0.227	0.179	0.179	0.178	0.185
0.4	0.344	0.344	0.36	0.353
0.25	0.223	0.233	0.233	0.234
0.234	0.186	0.192	0.193	0.206
0.438	0.459	0.433	0.469	0.449
0.306	0.306	0.296	0.315	0.297
0.342	0.329	0.313	0.342	0.337];

fmf = [0.466	0.458	0.472	0.5	0.42
0.292	0.29	0.301	0.31	0.273
0.36	0.353	0.365	0.396	0.32
0.474	0.424	0.445	0.483	0.45
0.328	0.294	0.301	0.325	0.31
0.364	0.291	0.306	0.338	0.337
0.54	0.552	0.547	0.569	0.554
0.376	0.374	0.367	0.4	0.366
0.469	0.482	0.47	0.48	0.482];

mmf = [ 0.452	0.44	0.433	0.465	0.383
0.277	0.282	0.285	0.303	0.243
0.332	0.331	0.315	0.357	0.29
0.452	0.422	0.413	0.435	0.439
0.306	0.281	0.294	0.306	0.274
0.336	0.274	0.284	0.301	0.308
0.521	0.557	0.538	0.509	0.517
0.368	0.369	0.354	0.364	0.347
0.453	0.462	0.446	0.436	0.458];

fcc = [0.582	0.633	0.569	0.468	0.576
0.353	0.345	0.338	0.328	0.333
0.548	0.555	0.516	0.509	0.499
0.564	0.585	0.557	0.603	0.599
0.351	0.376	0.345	0.359	0.365
0.511	0.51	0.502	0.531	0.532
0.571	0.564	0.549	0.587	0.567
0.379	0.363	0.35	0.366	0.36
0.509	0.489	0.472	0.51	0.507];

mcc = [0.449	0.619	0.528	0.557	0.524
0.325	0.342	0.322	0.33	0.321
0.519	0.545	0.494	0.497	0.483
0.578	0.57	0.542	0.573	0.555
0.344	0.364	0.342	0.35	0.349
0.501	0.493	0.482	0.501	0.497
0.549	0.545	0.52	0.55	0.525
0.365	0.35	0.343	0.349	0.341
0.486	0.464	0.442	0.488	0.474];

flame_stats.bp.low.acc    = fbp(1,:);
flame_stats.bp.low.f1     = fbp(2,:);
flame_stats.bp.low.auprc  = fbp(3,:);
flame_stats.bp.med.acc    = fbp(4,:);
flame_stats.bp.med.f1     = fbp(5,:);
flame_stats.bp.med.auprc  = fbp(6,:);
flame_stats.bp.high.acc   = fbp(7,:);
flame_stats.bp.high.f1    = fbp(8,:);
flame_stats.bp.high.auprc = fbp(9,:);
flame_stats.mf.low.acc    = fmf(1,:);
flame_stats.mf.low.f1     = fmf(2,:);
flame_stats.mf.low.auprc  = fmf(3,:);
flame_stats.mf.med.acc    = fmf(4,:);
flame_stats.mf.med.f1     = fmf(5,:);
flame_stats.mf.med.auprc  = fmf(6,:);
flame_stats.mf.high.acc   = fmf(7,:);
flame_stats.mf.high.f1    = fmf(8,:);
flame_stats.mf.high.auprc = fmf(9,:);
flame_stats.cc.low.acc    = fcc(1,:);
flame_stats.cc.low.f1     = fcc(2,:);
flame_stats.cc.low.auprc  = fcc(3,:);
flame_stats.cc.med.acc    = fcc(4,:);
flame_stats.cc.med.f1     = fcc(5,:);
flame_stats.cc.med.auprc  = fcc(6,:);
flame_stats.cc.high.acc   = fcc(7,:);
flame_stats.cc.high.f1    = fcc(8,:);
flame_stats.cc.high.auprc = fcc(9,:);

mash_stats.bp.low.acc    = mbp(1,:);
mash_stats.bp.low.f1     = mbp(2,:);
mash_stats.bp.low.auprc  = mbp(3,:);
mash_stats.bp.med.acc    = mbp(4,:);
mash_stats.bp.med.f1     = mbp(5,:);
mash_stats.bp.med.auprc  = mbp(6,:);
mash_stats.bp.high.acc   = mbp(7,:);
mash_stats.bp.high.f1    = mbp(8,:);
mash_stats.bp.high.auprc = mbp(9,:);
mash_stats.mf.low.acc    = mmf(1,:);
mash_stats.mf.low.f1     = mmf(2,:);
mash_stats.mf.low.auprc  = mmf(3,:);
mash_stats.mf.med.acc    = mmf(4,:);
mash_stats.mf.med.f1     = mmf(5,:);
mash_stats.mf.med.auprc  = mmf(6,:);
mash_stats.mf.high.acc   = mmf(7,:);
mash_stats.mf.high.f1    = mmf(8,:);
mash_stats.mf.high.auprc = mmf(9,:);
mash_stats.cc.low.acc    = mcc(1,:);
mash_stats.cc.low.f1     = mcc(2,:);
mash_stats.cc.low.auprc  = mcc(3,:);
mash_stats.cc.med.acc    = mcc(4,:);
mash_stats.cc.med.f1     = mcc(5,:);
mash_stats.cc.med.auprc  = mcc(6,:);
mash_stats.cc.high.acc   = mcc(7,:);
mash_stats.cc.high.f1    = mcc(8,:);
mash_stats.cc.high.auprc = mcc(9,:);
