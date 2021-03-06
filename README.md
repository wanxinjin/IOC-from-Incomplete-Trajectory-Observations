## Inverse Optimal Control from Incomplete Trajectory Observations

This project contains the source codes for the paper: Inverse Optimal Control from Incomplete Trajectory Observations (Published in the 
International Journal of Robotics Research) by Wanxin Jin, Dana Kulic, Shaoshuai Mou, and Sandra Hirche. Please find the paper at https://arxiv.org/abs/1803.07696.


### Folders
The current version of the project consists of three folders:

* **_OPTCON_** : contains multiple optimal control solvers to solve different optimal control problems.
* **_LQR_** : contains all source codes for  LQR examples in the paper. You can directly run each example.
* **_RobotArm_** : contains all source codes for  robot arm examples in the paper. You can directly run each example

Note that all the codes are implemented in MATLAB. Prior to using OPTCON to solve your optimal control problems, make sure you have installed CasADi in your MATLAB. 
Please refer to https://web.casadi.org/get/ for how to install CasADi in your MATLAB.




### How to Run the Codes
Directly run each script in LQR or RobotArm folder. Important lines in codes are commented. 


### Contact Information and Citation
If you have any question for the codes, please feel free to let me known via email:

   * name: wanxin jin (he/his)
   * email: wanxinjin@gmail.com


If you find this project helpful in your publications, please consider citing our paper.
    
    @article{jin2018inverse,
      title={Inverse optimal control with incomplete observations},
      author={Jin, Wanxin and Kuli{\'c}, Dana and Mou, Shaoshuai and Hirche, Sandra},
      journal={arXiv preprint arXiv:1803.07696},
      year={2018}
    }
 