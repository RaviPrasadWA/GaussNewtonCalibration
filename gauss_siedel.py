def calibrate_find_delta(dS, JS, delta):
    mu = 0.0
    for i in range(6):
        for j in range(i+1, 6):
            mu = JS[i][j]/float(JS[i][i])
            if( mu != 0.0 ):
                dS[j] -= mu*dS[i]
                for k in range(j,6):
                    JS[k][j] -= mu*JS[k][i]
                    
    for i in range(5,-1,-1):
        dS[i] /= float(JS[i][i])
        JS[i][i] = float(1.0)
        for j in range(i):
            mu = JS[i][j]
            dS[j] -= mu*dS[i]
            JS[i][j] = float(0.0)

    for i in range(6):
        delta[i] = dS[i]

def calibrate_update_matrices(dS, JS, beta, data):
    dx = 0
    b = 0
    residual = 1.0
    jacobian = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    
    for j in range(3):
        b = float(beta[3+j])
        dx = float(data[j]) - float(beta[j])
        residual -= b*b*dx*dx;
        jacobian[j] = float(2.0)*b*b*dx;
        jacobian[3+j] = float(-2.0)*b*dx*dx;

    for j in range(6):
        dS[j] += jacobian[j]*residual;
        for k in range(6):
            JS[j][k] += jacobian[j]*jacobian[k]
        

def calibrate_reset_matrices(dS, JS):
    for j in range(6):
        dS[j] = 0.0
        for k in range(6):
            JS[j][k] = 0.0

def accel_calib( accel_sample, accel_offsets, accel_scale):
    num_iterations  = 0
    eps             = 0.000000001
    change          = 100.0
    data            = [ 0.0, 0.0, 0.0]
    beta            = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    delta           = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ds              = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    JS              = []
    for i in range(6):
        JS.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    success         = True

    beta[0] = beta[1] = beta[2] = 0
    beta[3] = beta[4] = beta[5] = 1/float(9.78033)

    while( num_iterations < 20 and change > eps ):
        num_iterations += 1
        calibrate_reset_matrices(ds, JS)

        for i in range(6):
            data[0] = accel_sample[i][0] # 0 -> x
            data[1] = accel_sample[i][1] # 1 -> y
            data[2] = accel_sample[i][2] # 2 -> z
            calibrate_update_matrices(ds, JS, beta, data)
        calibrate_find_delta(ds, JS, delta)

        change = delta[0]*delta[0] +\
                 delta[0]*delta[0] +\
                 delta[1]*delta[1] +\
                 delta[2]*delta[2] +\
                 delta[3]*delta[3] / float(beta[3]*beta[3]) +\
                 delta[4]*delta[4] / float(beta[4]*beta[4]) +\
                 delta[5]*delta[5] / float(beta[5]*beta[5])


        for i in range(6) :
            beta[i] -= delta[i]
            
    accel_scale[0]          = round( beta[3] * float(9.78033), 2)
    accel_scale[1]          = round( beta[4] * float(9.78033), 2)
    accel_scale[2]          = round( beta[5] * float(9.78033), 2)
    accel_offsets[0]        = round( beta[0] * accel_scale[0], 2)
    accel_offsets[1]        = round( beta[1] * accel_scale[1], 2)
    accel_offsets[2]        = round( beta[2] * accel_scale[2], 2)

    valid_scale = False
    valid_offset = False

    if ( ( float(accel_scale[0]-1.0) > 0.1 ) or  (float(accel_scale[1]-1.0) > 0.1 ) or  (float(accel_scale[2]-1.0) > 0.1 )  ):
        valid_scale = True

    if ( ( float(accel_offsets[0]) > 3.5 ) or  (float(accel_offsets[1]) > 3.5 ) or  (float(accel_offsets[2]) > 3.5 )  ):
        valid_offset = True

    print valid_scale, accel_offsets

    if( valid_scale and valid_offset ):
        print '----------------------------'
        print 'Scaling -> ', accel_scale
        print 'Offsets -> ', accel_offsets
        print '----------------------------'


accel_sample_frame = [ [-0.11, 1.52, -9.04],
                       [0.22, -9.78, -0.45],
                       [0.79, 9.72, 0.80],
                       [9.84, 0.08, 0.95],
                       [-9.72, 0.20, -0.13],
                       [-0.16, -1.64, 10.47]
                       ]

accel_my_frame = [ [-0.96, 1.70, -8.98],
                   [0.05, -8.71, -2.08],
                   [0.35, 8.17, 1.27],
                   [9.09, 0.91, 0.36],
                   [-8.30, 0.60, -0.52],
                   [0.83, -1.56, 9.52] ]

accel_calib(accel_my_frame, [0,0,0], [0,0,0])
accel_calib(accel_sample_frame, [0,0,0], [0,0,0])
