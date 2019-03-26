#include <map_creator/kinematics.h>

#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;
#if IK_VERSION > 54
#define IKREAL_TYPE IkReal  // for IKFast 56,61
#else
#define IKREAL_TYPE IKReal  // for IKFast 54
#endif

IKREAL_TYPE eerot[9], eetrans[3];
namespace kinematics
{
float Kinematics::SIGN(float x)
{
    return (x >= 0.0f) ? +1.0f : -1.0f;
}

float Kinematics::NORM(float a, float b, float c, float d)
{
    return sqrt(a * a + b * b + c * c + d * d);
}

void Kinematics::getPoseFromFK(const std::vector< double > joint_values, std::vector< double >& pose)
{
#if IK_VERSION > 54
    // for IKFast 56,61
    unsigned int num_of_joints = GetNumJoints();
    unsigned int num_free_parameters = GetNumFreeParameters();
#else
    // for IKFast 54
    unsigned int num_of_joints = getNumJoints();
    unsigned int num_free_parameters = getNumFreeParameters();
#endif
    IKREAL_TYPE joints[num_of_joints];

    // cout<<joint_values[2]<<endl;

    for (unsigned int i = 0; i < num_of_joints; i++)
    {
        joints[i] = joint_values[i];
    }
#if IK_VERSION > 54
    // for IKFast 56,61
    ComputeFk(joints, eetrans, eerot);  // void return
#else
    // for IKFast 54
    fk(joints, eetrans, eerot);  // void return
#endif
    // cout<<"translation: "<<eetrans[0]<<" "<<eetrans[1]<<" "<<eetrans[2]<<endl;
    // Convert rotation matrix to quaternion (Daisuke Miyazaki)
    float q0 = (eerot[0] + eerot[4] + eerot[8] + 1.0f) / 4.0f;
    float q1 = (eerot[0] - eerot[4] - eerot[8] + 1.0f) / 4.0f;
    float q2 = (-eerot[0] + eerot[4] - eerot[8] + 1.0f) / 4.0f;
    float q3 = (-eerot[0] - eerot[4] + eerot[8] + 1.0f) / 4.0f;
    if (q0 < 0.0f)
        q0 = 0.0f;
    if (q1 < 0.0f)
        q1 = 0.0f;
    if (q2 < 0.0f)
        q2 = 0.0f;
    if (q3 < 0.0f)
        q3 = 0.0f;
    q0 = sqrt(q0);
    q1 = sqrt(q1);
    q2 = sqrt(q2);
    q3 = sqrt(q3);
    if (q0 >= q1 && q0 >= q2 && q0 >= q3)
    {
        q0 *= +1.0f;
        q1 *= SIGN(eerot[7] - eerot[5]);
        q2 *= SIGN(eerot[2] - eerot[6]);
        q3 *= SIGN(eerot[3] - eerot[1]);
    }
    else if (q1 >= q0 && q1 >= q2 && q1 >= q3)
    {
        q0 *= SIGN(eerot[7] - eerot[5]);
        q1 *= +1.0f;
        q2 *= SIGN(eerot[3] + eerot[1]);
        q3 *= SIGN(eerot[2] + eerot[6]);
    }
    else if (q2 >= q0 && q2 >= q1 && q2 >= q3)
    {
        q0 *= SIGN(eerot[2] - eerot[6]);
        q1 *= SIGN(eerot[3] + eerot[1]);
        q2 *= +1.0f;
        q3 *= SIGN(eerot[7] + eerot[5]);
    }
    else if (q3 >= q0 && q3 >= q1 && q3 >= q2)
    {
        q0 *= SIGN(eerot[3] - eerot[1]);
        q1 *= SIGN(eerot[6] + eerot[2]);
        q2 *= SIGN(eerot[7] + eerot[5]);
        q3 *= +1.0f;
    }
    else
    {
        printf("Error while converting to quaternion! \n");
    }
    float r = NORM(q0, q1, q2, q3);
    q0 /= r;
    q1 /= r;
    q2 /= r;
    q3 /= r;
    pose.push_back(eetrans[0]);
    pose.push_back(eetrans[1]);
    pose.push_back(eetrans[2]);
    pose.push_back(q1);
    pose.push_back(q2);
    pose.push_back(q3);
    pose.push_back(q0);
}

bool Kinematics::isIKSuccess(const std::vector< double >& pose, std::vector< double >& joints, int& numOfSolns, std::vector<bool>& configsFound)
{
#if IK_VERSION > 54
    // for IKFast 56,61
    unsigned int num_of_joints = GetNumJoints();
    unsigned int num_free_parameters = GetNumFreeParameters();
#else
    // for IKFast 54
    unsigned int num_of_joints = getNumJoints();
    unsigned int num_free_parameters = getNumFreeParameters();
#endif

#if IK_VERSION > 54
    // IKREAL_TYPE joints[num_of_joints];
    // for IKFast 56,61
    IkSolutionList< IKREAL_TYPE > solutions;
#else
    // for IKFast 54
    std::vector< IKSolution > vsolutions;
#endif
    std::vector< IKREAL_TYPE > vfree(num_free_parameters);
    eetrans[0] = pose[0];
    eetrans[1] = pose[1];
    eetrans[2] = pose[2];
    double qw = pose[6];
    double qx = pose[3];
    double qy = pose[4];
    double qz = pose[5];
    const double n = 1.0f / sqrt(qx * qx + qy * qy + qz * qz + qw * qw);
    qw *= n;
    qx *= n;
    qy *= n;
    qz *= n;
    eerot[0] = 1.0f - 2.0f * qy * qy - 2.0f * qz * qz; // T01
    eerot[1] = 2.0f * qx * qy - 2.0f * qz * qw;
    eerot[2] = 2.0f * qx * qz + 2.0f * qy * qw;
    eerot[3] = 2.0f * qx * qy + 2.0f * qz * qw;        // T11
    eerot[4] = 1.0f - 2.0f * qx * qx - 2.0f * qz * qz;
    eerot[5] = 2.0f * qy * qz - 2.0f * qx * qw;
    eerot[6] = 2.0f * qx * qz - 2.0f * qy * qw;        // T21
    eerot[7] = 2.0f * qy * qz + 2.0f * qx * qw;
    eerot[8] = 1.0f - 2.0f * qx * qx - 2.0f * qy * qy;
    // for(std::size_t i = 0; i < vfree.size(); ++i)
    // vfree[i] = atof(argv[13+i]);
    // TODO: the user have to define the number of free parameters for the manipulator if it has more than 6 joints. So
    // currently more than 6 joints are not supported yet.

#if IK_VERSION > 54 //(using this version > 54 :: vamsi)
    // for IKFast 56,61
    bool b1Success = ComputeIk(eetrans, eerot, vfree.size() > 0 ? &vfree[0] : NULL, solutions);

#else
    // for IKFast 54
    bool b2Success = ik(eetrans, eerot, vfree.size() > 0 ? &vfree[0] : NULL, vsolutions);

#endif

#if IK_VERSION > 54
    // for IKFast 56,61
    unsigned int num_of_solutions = (int)solutions.GetNumSolutions();
    numOfSolns = num_of_solutions;

#else
    // for IKFast 54
    unsigned int num_of_solutions = (int)vsolutions.size();
    numOfSolns = num_of_solutions;
#endif

    joints.resize(num_of_joints);

#if IK_VERSION > 54
    // for IKFast 56,61

    if (!b1Success)
    {
        return false;
    }
    else
    {
        //cout<<"Found ik solutions: "<< num_of_solutions<<endl;
        const IkSolutionBase< IKREAL_TYPE >& sol = solutions.GetSolution(0);
        int this_sol_free_params = (int)sol.GetFree().size();

        if( this_sol_free_params <= 0){
            sol.GetSolution(&joints[0], NULL);

            ///////change/////////////////////////////////////////////////////////////
            std::vector<IkReal> solvalues(GetNumJoints());
            for(std::size_t i = 0; i < solutions.GetNumSolutions(); ++i) {
                const IkSolutionBase<IkReal>& sol = solutions.GetSolution(i);
                //printf("sol%d (free=%d): ", (int)i, (int)sol.GetFree().size()); // printing the number of solution and free joints

                std::vector<IkReal> vsolfree(sol.GetFree().size());

                sol.GetSolution(&solvalues[0],vsolfree.size()>0?&vsolfree[0]:NULL);
                //                for( std::size_t j = 0; j < solvalues.size(); ++j)
                //                    printf("%.15f, ", solvalues[j]); //printing all the joint values of the solution
                //                printf("\n");
            }

            const double ZERO_THRESH = 0.00000001;
            const double PI = M_PI;
#define UR5_PARAMS
#ifdef UR5_PARAMS
            const double d1 =  0.089159;
            const double a2 = -0.42500;
            const double a3 = -0.39225;
            const double d4 =  0.10915;
            const double d5 =  0.09465;
            const double d6 =  0.0823;
#endif
            bool p1 = false, p3 = false, p5 =false;
            bool m1 = false, m3 = false, m5 =false;

            ////////// joint 1/////////////////
            double  q10, q11;
            double A = d6*eerot[4] - eerot[5];
            double B = d6*eerot[1] - eerot[2];
            double R = A*A + B*B;
            double arccos = acos(d4 / sqrt(R)) ;
            double arctan = atan2(-B, A);
            double pos = arccos + arctan;
            double neg = -arccos + arctan;
            if(fabs(pos) < ZERO_THRESH)
                pos = 0.0;
            if(fabs(neg) < ZERO_THRESH)
                neg = 0.0;
            if(pos >= 0.0)
                q10 = pos;
            else
                q10 = 2.0*PI + pos;
            if(neg >= 0.0)
                q11 = neg;
            else
                q11 = 2.0*PI + neg;

            if(std::fabs(solvalues[0]-q10) < 1e-6)
                p1 = true;    // right

            if(std::fabs(solvalues[0]-q11) < 1e-6)
                m1 = true;    // left

            /////// joint 5////////////
            double q500, q501, q510, q511;
            if (p1)
            {
                double numer = (eerot[2]*sin(q10) - eerot[5]*cos(q10)-d4);
                double div;
                if(fabs(fabs(numer) - fabs(d6)) < ZERO_THRESH)
                    div = SIGN(numer) * SIGN(d6);
                else
                    div = numer / d6;
                double arccos = acos(div);
                q500= arccos;
                q501 = 2.0*PI - arccos;
                if(std::fabs(solvalues[5]-q500) < 1e-6)
                    p1 = true;

                if(std::fabs(solvalues[5]-q501) < 1e-6)
                    m1 = true;
            }
            else if (m1)
            {
                double numer = (eerot[2]*sin(q11) - eerot[5]*cos(q11)-d4);
                double div;
                if(fabs(fabs(numer) - fabs(d6)) < ZERO_THRESH)
                    div = SIGN(numer) * SIGN(d6);
                else
                    div = numer / d6;
                double arccos = acos(div);
                q510= arccos;
                q511 = 2.0*PI - arccos;
                if(std::fabs(solvalues[5]-q510) < 1e-6)
                    p5 = true;

                if(std::fabs(solvalues[5]-q511) < 1e-6)
                    m5 = true;
            }
            /////// joint 3////////////
            double q3 = fabs(solvalues[3]);
            if(solvalues[3] == q3)
                p3 = true;  // elbow up

            if(solvalues[3] < q3)
                m3 = true;

            if (p1 == true && p3 == true && p5== true)
                configsFound[0] = true;
            if (p1 == true && m3 == true && p5== true)
                configsFound[1] = true;
            if (p1 == true && m3 == true && m5== true)
                configsFound[2] = true;
            if (p1 == true && p3 == true && m5== true)
                configsFound[3] = true;
            if (m1 == true && p3 == true && p5== true)
                configsFound[4] = true;
            if (m1 == true && m3 == true && p5== true)
                configsFound[5] = true;
            if (m1 == true && m3 == true && m5== true)
                configsFound[6] = true;
            if (m1 == true && p3 == true && m5== true)
                configsFound[7] = true;
            ////////////////////////////////////////////////////////////////
        }
        else{
            static std::vector< IKREAL_TYPE > vsolfree;
            vsolfree.resize(this_sol_free_params);
            sol.GetSolution(&joints[0], &vsolfree[0]);
        }
        return true;
    }
#else
    if (!b2success)
    {
        return false;
    }
    else
    {
        cout<<"Found ik solutions: "<< num_of_solutions<<endl;
        int this_sol_free_params = (int)vsolutions[i].GetFree().size();
        std::vector< IKREAL_TYPE > vsolfree(this_sol_free_params);
        vsolutions[i].GetSolution(&solvalues[0], vsolfree.size() > 0 ? &vsolfree[0] : NULL);
        return true;
    }
#endif
}

const std::string Kinematics::getRobotName()
{
    const char* hash = GetKinematicsHash();

    std::string part = hash;
    part.erase(0, 22);
    std::string name = part.substr(0, part.find(" "));
    return name;
}

bool Kinematics::isIkSuccesswithTransformedBase(const geometry_msgs::Pose& base_pose,
                                                const geometry_msgs::Pose& grasp_pose, std::vector<double>& joint_soln,int& numOfSolns)
{
    // Creating a transformation out of base pose
    tf2::Vector3 base_vec(base_pose.position.x, base_pose.position.y, base_pose.position.z);
    tf2::Quaternion base_quat(base_pose.orientation.x, base_pose.orientation.y, base_pose.orientation.z,
                              base_pose.orientation.w);
    base_quat.normalize();
    tf2::Transform base_trns;
    base_trns.setOrigin(base_vec);
    base_trns.setRotation(base_quat);

    // Inverse of the transformation
    tf2::Transform base_trns_inv;
    base_trns_inv = base_trns.inverse();

    // Creating a transformation of grasp pose
    tf2::Vector3 grasp_vec(grasp_pose.position.x, grasp_pose.position.y, grasp_pose.position.z);
    tf2::Quaternion grasp_quat(grasp_pose.orientation.x, grasp_pose.orientation.y, grasp_pose.orientation.z,
                               grasp_pose.orientation.w);
    grasp_quat.normalize();
    tf2::Transform grasp_trns;
    grasp_trns.setOrigin(grasp_vec);
    grasp_trns.setRotation(grasp_quat);

    // Transforming grasp pose to origin from where we can check for Ik
    tf2::Transform new_grasp_trns;
    // new_grasp_trns = grasp_trns * base_trns_inv;
    new_grasp_trns = base_trns_inv * grasp_trns;
    // Creating a new grasp pose in the origin co-ordinate
    std::vector< double > new_grasp_pos;
    tf2::Vector3 new_grasp_vec;
    tf2::Quaternion new_grasp_quat;
    new_grasp_vec = new_grasp_trns.getOrigin();
    new_grasp_quat = new_grasp_trns.getRotation();
    new_grasp_quat.normalize();
    new_grasp_pos.push_back(new_grasp_vec[0]);
    new_grasp_pos.push_back(new_grasp_vec[1]);
    new_grasp_pos.push_back(new_grasp_vec[2]);
    new_grasp_pos.push_back(new_grasp_quat[0]);
    new_grasp_pos.push_back(new_grasp_quat[1]);
    new_grasp_pos.push_back(new_grasp_quat[2]);
    new_grasp_pos.push_back(new_grasp_quat[3]);

    // Check the new grasp_pose for Ik
    Kinematics k;
    //std::vector< double > joints;

    //joints.resize(6);
    if (k.isIKSuccess(new_grasp_pos,  joint_soln, numOfSolns))
        return true;
    else
        return false;
}





};
