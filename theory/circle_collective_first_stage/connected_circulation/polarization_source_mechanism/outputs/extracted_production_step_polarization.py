def archived_step_polarization(current, para, noise_increment, dis_last):
    state = np.empty((3, para['num'], 2))
    state[:, :, 0] = current
    rDotList = np.empty((2, para['num'], 2))
    count = 0
    CurState = state[:, :, count]
    CurLocation = state[0:2, :, count]
    CurHeading = state[2, :, count]
    ex = np.cos(CurHeading)
    ey = np.sin(CurHeading)
    e = np.vstack([ex, ey])
    [U, wHydro] = GetHydro(CurState, para)
    [U_BEM, w_BEM] = GetHydroBEM(CurState, para)
    voro = GetVoronoi(CurLocation, para)
    wVision = GetVision(voro, CurState, para)
    wNoise = noise_increment
    wWall = GetWallAvoid(CurState, para)
    thetaDot = wHydro + w_BEM + para['light'] * (wWall + wVision)
    thetaDot = np.clip(thetaDot, -np.pi / para['dt'], np.pi / para['dt'])
    rDot = e + U + U_BEM
    if_jump = Check_jumping(CurState, rDot, para)
    thetaDot_nohydro = para['light'] * (wWall + wVision)
    thetaDot_nohydro = np.clip(thetaDot_nohydro, -np.pi / para['dt'], np.pi / para['dt'])
    if if_jump.any():
        rDot[:, if_jump] = e[:, if_jump]
        thetaDot[:, if_jump] = thetaDot_nohydro[:, if_jump]
    if np.any(abs(rDot) > 500):
        where_exceed = np.unique(np.where(abs(rDot) > 500)[1])
        rDot[:, where_exceed] = e[:, where_exceed]
        thetaDot[:, where_exceed] = thetaDot_nohydro[:, where_exceed]
    count += 1
    state[:, :, count] = CurState + np.vstack([rDot, thetaDot]) * para['dt'] + wNoise
    dis = np.zeros(para['num'])
    for i in range(para['num']):
        dis[i] = cal_signed_distance_function(state[0:2, i, count])
    if any(dis - safedis >= 0):
        locationCheck = np.sign(dis - safedis)
        if dis_last is None:
            warningFish = np.where(locationCheck >= 0)[0]
        else:
            warningFish = np.where(np.logical_and(locationCheck >= 0, dis_last < dis))[0]
        for k in range(len(warningFish)):
            rDot[:, warningFish[k]] = np.zeros(2)
        state[:, :, count] = CurState + np.vstack([rDot, thetaDot]) * para['dt'] + wNoise
    dis_last = np.copy(dis)
    rDotList[:, :, count] = rDot
    angular = np.vstack([np.asarray(wHydro).ravel(), np.asarray(w_BEM).ravel(), np.asarray(wWall).ravel(), np.asarray(wVision).ravel(), np.asarray(thetaDot).ravel()])
    flags = dict(wall_clip=int(np.sum(abs(wWall) >= np.pi / para['dt'])), angular_clip=int(np.sum(abs(thetaDot) >= np.pi / para['dt'])), collision=int(if_jump.sum()), safety_candidates=int(np.sum(dis >= safedis)))
    return (state[:, :, count].copy(), rDot.copy(), dis_last.copy(), flags, angular, voro.copy())
