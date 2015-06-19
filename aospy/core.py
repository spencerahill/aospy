# def _aospy_inst(proj=False, model=False, run=False, var=False):
#     """Convert string matching aospy object names to class instances."""
#     def _to_list(obj):
#         if type(obj) is str:
#             list_obj = [obj]
#         else:
#             list_obj = obj
#         return list_obj

#     pr, md, rn, vr = [], [], [], []
#     for p, m, r, v in zip(_to_list(proj), _to_list(model), _to_list(run),
#                           _to_list(var)):
#         if p:
#             p = _proj_inst(p)
#             pr.append(p)
#         if m:
#             m = _model_inst(m, p)
#             md.append(m)
#         if r:
#             r = _run_inst(r, m, p)
#             rn.append(r)
#         if v:
#             vr.append(_var_inst(v))

#     def _strip_singleton_list(l):
#         """Strip extra list layers."""
#         while True:
#             try:
#                 l = l[0]
#             except TypeError:
#                 break
#         return l

#     if type(proj) is str or len(proj) == 1:
#         pr = _strip_singleton_list(pr)
#     if type(model) is str or len(model) == 1:
#         md = _strip_singleton_list(md)
#     if type(run) is str or len(run) == 1:
#         rn = _strip_singleton_list(rn)
#     if type(var) is str or len(var) == 1:
#         vr = _strip_singleton_list(vr)

#     return pr, md, rn, vr
