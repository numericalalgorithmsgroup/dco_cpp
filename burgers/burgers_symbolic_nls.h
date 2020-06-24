#pragma once

void newton_adjoint(DCO_EAO_T* D) {
  const vector<DCO_BT>& y=D->read_data<vector<DCO_BT>>();
  vector<DCO_BT> ya(y.size());
  D->get_output_adjoint(ya);
  vector<DCO_BT> A((y.size()-2)*3+4,0);
  dfdy(y,A,/*transpose=*/true);
  LU(A); FS(A,ya); BS(A,ya);
  D->increment_input_adjoint(ya);
}

void newton(const vector<DCO_T>& yp, vector<DCO_T>& y) {
  DCO_TAPE_T* tape=dco::tape(yp);
  DCO_EAO_T* D=tape->create_callback_object<DCO_EAO_T>();
  vector<DCO_BT> ypv=D->register_input(yp);
  vector<DCO_BT> yv=dco::value(y);
  newton(ypv,yv);
  D->write_data(yv);
  y=D->register_output(yv); 
  tape->insert_callback(newton_adjoint, D);
}
  

