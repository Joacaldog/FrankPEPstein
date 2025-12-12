
import py3Dmol
v = py3Dmol.view()
print(dir(v))
try:
    print("HAS _make_html:", hasattr(v, '_make_html'))
except:
    pass
