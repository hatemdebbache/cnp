Attribute VB_Name = "DERIVE"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << DERIVE >>
'  Requis   :           << EVAL >>
'  Externes :              Fonc(nf,x)
'
'  Descrip. : Dérivation numérique
'
Option Explicit
Sub Tchebychev1(ByVal ordre As Integer, _
                ByVal M As Integer, _
                ByRef C() As Double, _
                ByRef d() As Double)

    ' c(i-1) = c(i+1) +2*(i-1)*c(i)
End Sub
Sub Tchebychev2(ByVal ordre As Integer, _
                ByVal M As Integer, _
                ByRef C() As Double)

    ' c(i-1) = c(i+1) +2*(i-1)*c(i)
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'                   Dérivées de Polynômes
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Function Horner1(ByVal Deg As Integer, _
                 ByRef Coef() As Double, _
                 ByVal x As Double, _
        Optional ByRef Pol As Double, _
        Optional ByVal BASE As Integer = 1, _
        Optional ByVal SENS As Integer = 1) As Double

  ' Evaluation de P'(x) par un schéma de Horner
  '
  '
  '
  
  ' Déclaration
    Dim k As Integer, i As Integer
    Dim Der As Double
    
  ' Contrôle
    If Deg < 0 Then Stop         ' Données non valides
    If Abs(SENS) <> 1 Then Stop  ' Données non valides
    
  ' Début du programme
    If SENS = 1 Then k = Deg + BASE Else k = BASE
    Der = 0
    Pol = Coef(k)
    For i = 1 To Deg
        k = k - SENS
        Der = x * Der + Pol
        Pol = x * Pol + Coef(k)
    Next i
    Horner1 = Der
End Function
Function Horner2(ByVal Deg As Integer, _
                 ByRef Coef() As Double, _
                 ByVal x As Double, _
        Optional ByRef Der As Double, _
        Optional ByRef Pol As Double, _
        Optional ByVal BASE As Integer = 1, _
        Optional ByVal SENS As Integer = 1) As Double

  ' Evaluation de P"(x) par un schéma de Horner
  '
  '
  '
  
  ' Déclaration
    Dim k As Integer, i As Integer
    Dim Sec As Double
    
  ' Contrôle
    If Deg < 0 Then Stop         ' Données non valides
    If Abs(SENS) <> 1 Then Stop  ' Données non valides
    
  ' Début du programme
    If SENS = 1 Then k = Deg + BASE Else k = BASE
    Sec = 0#
    Der = 0#
    Pol = Coef(k)
    For i = 1 To Deg
        k = k - SENS
        Sec = x * Sec + Der
        Der = x * Der + Pol
        Pol = x * Pol + Coef(k)
    Next i
    Horner2 = 2# * Sec
End Function
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'         Dérivées d'un polynôme d'interpolation
'             mis sous la forme de Newton
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Function Newton1(ByVal k1 As Integer, _
                 ByVal k2 As Integer, _
                 ByRef x() As Double, _
                 ByRef C() As Double, _
                 ByVal xo As Double, _
                 ByRef Der As Double) As Double

    Dim k As Integer, ko As Integer, N As Integer
    Dim Pol As Double, u As Double
    
    N = k2 - k1 + 1
    If (N < 1) Then Stop ' Données non valides
    Der = 0#
    If N = 1 Then Newton1 = C(1): Exit Function
    Pol = C(N)
    ko = k2 - N
    For k = N - 1 To 1 Step -1
        u = xo - x(k + ko)
        Der = u * Der + Pol
        Pol = u * Pol + C(k)
    Next k
    Newton1 = Pol
End Function
Function Newton2(ByVal k1 As Integer, _
                 ByVal k2 As Integer, _
                 ByRef x() As Double, _
                 ByRef C() As Double, _
                 ByVal xo As Double, _
                 ByRef Der As Double, _
                 ByRef Sec As Double) As Double

    Dim k As Integer, ko As Integer, N As Integer
    Dim Pol As Double, u As Double
    
    N = k2 - k1 + 1
    If (N < 1) Then Stop ' Données non valides
    Der = 0#: Sec = 0#
    If N = 1 Then Newton2 = C(1): Exit Function
    Pol = C(N)
    ko = k2 - N
    For k = N - 1 To 1 Step -1
        u = xo - x(k + ko)
        Sec = u * Sec + Der
        Der = u * Der + Pol
        Pol = u * Pol + C(k)
    Next k
    Sec = 2# * Sec
    Newton2 = Pol
End Function
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'                Formules de différences finies
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Function FDF(ByVal nf As Integer, _
             ByVal x As Double, _
             ByVal h As Double) As Double
             
    FDF = (Fonc(nf, x + h) - Fonc(nf, x)) / h
End Function
Function FDFC(ByVal L As Integer, _
              ByVal nf As Integer, _
              ByVal x As Double, _
              ByVal h As Double) As Double
             
    Select Case L
    Case 1  ' Dérivée première
        FDFC = 0.5 * (Fonc(nf, x + h) - Fonc(nf, x - h)) / h
    Case 2  ' Dérivée seconde
        FDFC = (Fonc(nf, x + h) - 2 * Fonc(nf, x) + Fonc(nf, x - h)) / (h * h)
    Case Else
        Stop
    End Select
End Function
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'             Méthodes itératives (Richardson)
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Function Richardson0(ByVal L As Integer, _
                     ByVal nf As Integer, _
                     ByVal x As Double, _
                     ByVal ho As Double, _
                     ByVal eps As Double, _
                     ByRef Der As Double, _
            Optional ByRef dif As Double) As Integer

    Const kmax As Integer = 20
    Dim k As Integer
    Dim h As Double, Ader As Double, Adif As Double
    
    h = ho
    Der = DERIVE.FDFC(L, nf, x, h)
    dif = 1E+30
    k = 0
    Do
        k = k + 1
        If k > kmax Then Richardson0 = 0: Exit Function
        Ader = Der
        Adif = dif
        h = 0.5 * h
        Der = DERIVE.FDFC(L, nf, x, h)
        dif = Abs(Der - Ader)
        If dif < eps Then Richardson0 = k: Exit Function
    Loop Until 3 * dif > Adif
    Der = Ader
    dif = Adif
    Richardson0 = -k
End Function
Sub DERIVE_Richardson1(ByVal MDF As Integer, _
                       ByVal nf As Integer, _
                       ByVal h As Double, _
                       ByVal x As Double, _
                       ByVal eps As Double, _
                       ByRef Der As Double, _
                       ByRef err As Integer)

    Const kmax As Integer = 15
    
    Dim d As Double, ancien As Double
    Dim k As Integer
    
    d = DERIVE_DFcentrale(MDF, nf, h, x)
    k = 0
    Do
        k = k + 1
        If (k > kmax) Then err = mNONCNV: Exit Sub
        h = 0.5 * h
        ancien = d
        d = DERIVE_DFcentrale(MDF, nf, h, x)
        Der = (4# * d - ancien) / 3#
    Loop Until Abs(Der - ancien) < eps
    err = k
End Sub
Sub DERIVE_Richardson2(ByVal MDF As Integer, _
                       ByVal nf As Integer, _
                       ByVal h As Double, _
                       ByVal x As Double, _
                       ByVal eps As Double, _
                       ByRef Der As Double, _
                       ByRef err As Integer)

    Const kmax As Integer = 15
    
    Dim r(0 To kmax) As Double
    Dim d As Double, ancien As Double, p As Double
    Dim k As Integer, j As Integer
    
    ancien = 1E+30
    r(0) = DERIVE_DFcentrale(MDF, nf, h, x)
    k = 0
    Do
        k = k + 1
        If (k > kmax) Then err = mNONCNV: Exit Sub
        h = 0.5 * h
        d = r(0)
        r(0) = DERIVE_DFcentrale(MDF, nf, h, x)
        p = 1
        For j = 1 To k
            p = 4# * p
            ancien = d
            d = r(j)
            r(j) = (p * r(j - 1) - ancien) / (p - 1#)
        Next j
    Loop Until Abs(r(k) - ancien) < eps
    Der = r(k)
    err = k
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'          Dérivées d'une fonction Spline cubique
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Function Spline1(ByVal N As Integer, _
                 ByRef x() As Double, _
                 ByRef y() As Double, _
                 ByRef s() As Double, _
                 ByVal xo As Double, _
                 ByRef Der As Double) As Double

    Dim k1 As Integer, k2 As Integer
    Dim v As Double, w As Double
    Dim v2 As Double, w2 As Double
    Dim u1 As Double, u2 As Double
    Dim h As Double, h6 As Double
    
    Call COEFF.Morceau(N, x(), xo, k1, k2, h, v, w, u1, u2)
    h6 = h / 6: v2 = v * v: w2 = w * w
    
    Der = (y(k2) - y(k1)) / h + _
          ((3# * v2 - 1) * s(k2) - _
           (3# * w2 - 1) * s(k1)) * h6
    
    Spline1 = v * y(k2) + w * y(k1) + _
             ((v2 - 1#) * u1 * s(k2) + _
              (w2 - 1#) * u2 * s(k1)) * h6
End Function
Function Spline2(ByVal N As Integer, _
                 ByRef x() As Double, _
                 ByRef y() As Double, _
                 ByRef s() As Double, _
                 ByVal xo As Double, _
                 ByRef Der As Double, _
                 ByRef Sec As Double) As Double

    Dim k1 As Integer, k2 As Integer
    Dim v As Double, w As Double
    Dim v2 As Double, w2 As Double
    Dim u1 As Double, u2 As Double
    Dim h As Double, h6 As Double
    
    Call COEFF.Morceau(N, x(), xo, k1, k2, h, v, w, u1, u2)
    h6 = h / 6: v2 = v * v: w2 = w * w
    
    Sec = v * s(k2) + w * s(k1)
    
    Der = (y(k2) - y(k1)) / h + _
         ((3# * v2 - 1) * s(k2) - _
          (3# * w2 - 1) * s(k1)) * h6
    
    Spline2 = v * y(k2) + w * y(k1) + _
             ((v2 - 1#) * u1 * s(k2) + _
              (w2 - 1#) * u2 * s(k1)) * h6
 End Function

