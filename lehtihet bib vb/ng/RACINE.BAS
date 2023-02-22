Attribute VB_Name = "RACINE"
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
'
'                     Module RACINE
'
' Requis :
'
'         fonction   Fonc(nf,x,[d],[s])
'
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
Private Const mKMAX As Integer = 30
Private Const mINFINI As Double = 1E+30
Private Const mNONVAL As Integer = -1
Private Const mNONCNV As Integer = -2
Private Const mSINGUL As Integer = -3
Private Const mDIVERG As Integer = -4
Public TabX(0 To mKMAX) As Double
Public TabY(0 To mKMAX) As Double
Public TabDelta(0 To mKMAX) As Double
Public TabLamda(0 To mKMAX) As Double
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
'
'                   Méthodes exactes
'
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
Public Sub RACINE_ESD(ByVal b As Double, _
                      ByVal c As Double, _
                      ByRef x1 As Double, _
                      ByRef x2 As Double, _
                      ByRef Nsol As Integer)

    ' Résolution d'une équation du second de degré avec
    ' réduction des effets dûs aux phénomènes de compensation.
    '
    ' A partir des coefficients b, c, la procédure détermine
    ' Nsol le nombre de racines réelles de x² + bx + c = 0.
    ' x1 et x2 sont les solutions recherchées.
    '
    ' Si Nsol = 2, x1 sera la solution ayant la plus petite
    '              valeur absolue.
    ' Si Nsol = 1, x1 sera la racine double
    ' Si Nsol = 0, x1 sera la partie réelle et x2 la partie
    '              imaginaire des solutions complexes conjuguées.
    
    Dim delta As Double
    
    x1 = -0.5 * b
    delta = x1 * x1 - c
    Nsol = 1 + Sgn(delta)
    x2 = Sqr(Abs(delta))
    If (delta > 0) Then
        If (x1 <> 0) Then x2 = x1 + Sgn(x1) * x2
        x1 = c / x2
    End If
 End Sub
Public Sub RACINE_ESD0(ByVal b As Double, _
                       ByVal c As Double, _
                       ByRef x1 As Double, _
                       ByRef x2 As Double, _
                       ByRef err As Integer)

    Dim delta As Double
    
    x2 = 0#
    x1 = -0.5 * b
    delta = x1 * x1 - c
    If delta = 0 Then
        err = 1
    ElseIf delta < 0 Then
        x2 = Sqr(-delta)
        err = 0
    Else
        x2 = Sqr(delta)
        If x1 > 0 Then x2 = x1 + x2 Else x2 = x1 - x2
        x1 = c / x2
        err = 2
    End If
 End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
'
'               Méthodes itératives fermées
'
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
Sub RACINE_Resultats(ByVal err As Integer)

    Dim i As Integer, U As Double
    
    For i = 0 To err
        U = Abs(TabX(i) - 0.5 * 3.1415926)
        Debug.Print i, Format(U, "0.000e+00")
  '      TabY(i) = Log(u + 1E-30) / Log(10)
  '      TabX(i) = i
    Next i
    Debug.Print
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
'
'               Méthodes itératives fermées
'
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
Public Sub RACINE_FER_Bolzano(ByVal nf As Integer, _
                              ByVal a As Double, _
                              ByVal b As Double, _
                              ByVal dx As Double, _
                              ByRef X As Double, _
                              ByRef err As Integer)
                                 
    Dim k As Integer
    Dim Y As Double, L As Double
    Dim xn As Double, xp As Double
    Dim yn As Double, yp As Double
    
    If prv_Cherche(nf, a, b, xn, xp, yn, yp) Then err = mNONVAL: Exit Sub
    L = xp - xn
    k = -1
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        L = 0.5 * L
        X = xp - L: TabX(k) = X
        Y = Fonc(nf, X)
        If (Y > 0#) Then xp = X
    Loop Until ((Abs(L) < dx) Or (Y = 0))
    err = k + 3
End Sub
Public Sub RACINE_FER_RegulaFalsi(ByVal nf As Integer, _
                                  ByVal a As Double, _
                                  ByVal b As Double, _
                                  ByVal dx As Double, _
                                  ByVal dy As Double, _
                                  ByRef X As Double, _
                                  ByRef err As Integer)
                             
    Dim k As Integer
    Dim xn As Double, xp As Double
    Dim yn As Double, yp As Double
    Dim Y As Double, xx As Double
    
    If prv_Cherche(nf, a, b, xn, xp, yn, yp) Then err = mNONVAL: Exit Sub
    X = xp
    k = -1
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        xx = X
        If (yp + yn) < 0 Then
            X = xp - (xp - xn) * (yp / (yp - yn)): TabX(k) = X
        Else
            X = xn - (xp - xn) * (yn / (yp - yn)): TabX(k) = X
        End If
        Y = Fonc(nf, X)
        If (Y > 0#) Then
            xp = X
            yp = Y
        Else
            xn = X
            yn = Y
        End If
        TabLamda(k) = xp - xn
        TabDelta(k) = X - xx
    Loop Until prv_Arret(dx, dy, X - xx, Y)
    err = k + 3
End Sub
Public Sub RACINE_FER_Hamming1(ByVal nf As Integer, _
                               ByVal a As Double, _
                               ByVal b As Double, _
                               ByVal dx As Double, _
                               ByVal dy As Double, _
                               ByRef X As Double, _
                               ByRef err As Integer)
                           
    Dim k As Integer
    Dim xn As Double, xp As Double
    Dim yn As Double, yp As Double
    Dim Y As Double, xx As Double
    
    If prv_Cherche(nf, a, b, xn, xp, yn, yp) Then err = mNONVAL: Exit Sub
    X = xp
    k = -1
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        xx = X
        X = xp - (xp - xn) * (yp / (yp - yn)): TabX(k) = X
        Y = Fonc(nf, X)
        If (Y > 0#) Then
            yn = 0.5 * yn
            xp = X
            yp = Y
        Else
            yp = 0.5 * yp
            xn = X
            yn = Y
        End If
        TabLamda(k) = xp - xn
        TabDelta(k) = X - xx
    Loop Until prv_Arret(dx, dy, X - xx, Y)
    err = k + 3
End Sub
Public Sub RACINE_FER_Hamming2(ByVal nf As Integer, _
                               ByVal a As Double, _
                               ByVal b As Double, _
                               ByVal dx As Double, _
                               ByVal dy As Double, _
                               ByRef X As Double, _
                               ByRef err As Integer)
                           
    Dim k As Integer, s As Integer
    Dim xn As Double, xp As Double
    Dim yn As Double, yp As Double
    Dim Y As Double, xx As Double
    
    If prv_Cherche(nf, a, b, xn, xp, yn, yp) Then err = mNONVAL: Exit Sub
    X = xp
    s = 0
    k = -1
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        xx = X
        X = xp - (xp - xn) * (yp / (yp - yn)): TabX(k) = X
        Y = Fonc(nf, X)
        If (Y > 0#) Then
            If (s > 0) Then yn = 0.5 * yn
            s = 1
            xp = X
            yp = Y
        Else
            If (s < 0) Then yp = 0.5 * yp
            s = -1
            xn = X
            yn = Y
        End If
        TabLamda(k) = xp - xn
        TabDelta(k) = X - xx
    Loop Until prv_Arret(dx, dy, X - xx, Y)
    err = k + 3
End Sub
Public Sub RACINE_FER_Hamming3(ByVal nf As Integer, _
                               ByVal a As Double, _
                               ByVal b As Double, _
                               ByVal dx As Double, _
                               ByVal dy As Double, _
                               ByRef X As Double, _
                               ByRef err As Integer)
    
    Dim k As Integer, s As Integer
    Dim xn As Double, xp As Double
    Dim yn As Double, yp As Double
    Dim Y As Double, xx As Double
    
    If prv_Cherche(nf, a, b, xn, xp, yn, yp) Then err = mNONVAL: Exit Sub
    X = xp
    s = 0
    k = -1
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        xx = X
        X = xp - (xp - xn) * (yp / (yp - yn)): TabX(k) = X
        Y = Fonc(nf, X)
        If (Y > 0#) Then
            If (s > 0) Then yn = yn * (yp / (yp + Y))
            s = 1
            xp = X
            yp = Y
        Else
            If (s < 0) Then yp = yp * (yn / (yn + Y))
            s = -1
            xn = X
            yn = Y
        End If
        TabLamda(k) = xp - xn
        TabDelta(k) = X - xx
    Loop Until prv_Arret(dx, dy, X - xx, Y)
    err = k + 3
End Sub
Public Sub RACINE_FER_Ridders(ByVal nf As Integer, _
                              ByVal a As Double, _
                              ByVal b As Double, _
                              ByVal dx As Double, _
                              ByVal dy As Double, _
                              ByRef X As Double, _
                              ByRef err As Integer)
                             
    Dim k As Integer
    Dim xn As Double, xp As Double
    Dim yn As Double, yp As Double
    Dim Y As Double, xx As Double, d As Double
    
    If prv_Cherche(nf, a, b, xn, xp, yn, yp) Then err = mNONVAL: Exit Sub
    X = 0.5 * (xp + xn)
    k = -1
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        xx = X
        X = 0.5 * (xp + xn): TabX(k) = X
        Y = Fonc(nf, X)
        d = Y * Y - yp * yn
        If (Y = 0) Then err = 2 * k + 3: Exit Sub
        If (Y > 0) Then
            xp = X
            yp = Y
        Else
            xn = X
            yn = Y
        End If
        X = X - (xp - xn) * (Y / Sqr(d)): TabX(k) = X
        Y = Fonc(nf, X)
        If (Y > 0) Then
            xp = X
            yp = Y
        Else
            xn = X
            yn = Y
        End If
        TabLamda(k) = xp - xn
        TabDelta(k) = X - xx
    Loop Until prv_Arret(dx, dy, X - xx, Y)
    err = 2 * k + 4
End Sub
Public Sub RACINE_FER_Muller(ByVal nf As Integer, _
                             ByVal a As Double, _
                             ByVal b As Double, _
                             ByVal dx As Double, _
                             ByVal dy As Double, _
                             ByRef X As Double, _
                             ByRef err As Integer)
                             
    Dim k As Integer, s As Integer
    Dim Y As Double, delta As Double
    Dim xn As Double, xp As Double
    Dim yn As Double, yp As Double
    Dim hn As Double, hp As Double
    Dim un As Double, up As Double
    Dim ALPHA As Double, BETA As Double
    Dim GAMMA As Double, w As Double
    
    If prv_Cherche(nf, a, b, xn, xp, yn, yp) Then err = mNONVAL: Exit Sub
    X = 0.5 * (xp + xn)
    Y = Fonc(nf, X)
    s = 0
    k = -1
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        un = xn - X: hn = (yn - Y) / un
        up = xp - X: hp = (yp - Y) / up
        If Y > 0 Then
            If (s > 0) Then yn = yn * (yp / (yp + Y))
            s = 1
            xp = X
            yp = Y
        Else
            If (s < 0) Then yp = yp * (yn / (yn + Y))
            s = -1
            xn = X
            yn = Y
        End If
        ALPHA = hp - hn
        BETA = 0.5 * (up * hn - un * hp)
        GAMMA = Y * (up - un)
        w = Sqr(BETA * BETA - ALPHA * GAMMA)
        If BETA < 0 Then w = -w
        w = BETA + w
        delta = GAMMA / w
        If (up + delta) * (un + delta) > 0 Then delta = w / ALPHA
        X = X - delta: TabX(k) = X
        Y = Fonc(nf, X)
        TabLamda(k) = xp - xn
        TabDelta(k) = delta
    Loop Until prv_Arret(dx, dy, delta, Y)
    err = k + 4
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
'
'               Méthodes itératives ouvertes
'
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
Public Sub RACINE_OUV_Newton(ByVal nf As Integer, _
                             ByVal a As Double, _
                             ByVal b As Double, _
                             ByVal dx As Double, _
                             ByVal dy As Double, _
                             ByVal xo As Double, _
                             ByRef X As Double, _
                             ByRef err As Integer, _
                    Optional ByVal M As Integer = 1)
                             
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim d As Double
Dim c As String
    
frm.List.Clear
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    Y = Fonc(nf, X, d)
    k = 0: TabX(k) = X
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        If (d = 0#) Then err = mSINGUL: Exit Sub
        delta = M * (Y / d)
c = Format(X, "0.000000    ") + Format(delta, "0.000000 ")
frm.List.AddItem (c)
        X = X - delta: TabX(k) = X
        If prv_Divergence(a, X, b) Then err = mDIVERG: Exit Sub
        Y = Fonc(nf, X, d)
        TabLamda(k) = Y
        TabDelta(k) = delta
    Loop Until prv_Arret(dx, dy, delta, Y)
    err = k + 1
c = Format(X, "0.000000 ")
frm.List.AddItem (c)
End Sub
Public Sub RACINE_OUV_VonMises(ByVal nf As Integer, _
                               ByVal a As Double, _
                               ByVal b As Double, _
                               ByVal dx As Double, _
                               ByVal dy As Double, _
                               ByVal xo As Double, _
                               ByVal h As Double, _
                               ByRef X As Double, _
                               ByRef err As Integer, _
                      Optional ByVal M As Integer = 1)
                            
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim d As Double
    
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    Y = Fonc(nf, X, d)
  ' d = Fonc(nf, x + h) - y
    If (d = 0#) Then err = mSINGUL: Exit Sub
    d = M / d
    k = 0: TabX(k) = X
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        delta = Y * d
        X = X - delta: TabX(k) = X
        If prv_Divergence(a, X, b) Then err = mDIVERG: Exit Sub
        Y = Fonc(nf, X)
        TabLamda(k) = Y
        TabDelta(k) = delta
    Loop Until prv_Arret(dx, dy, delta, Y)
    err = k + 1
End Sub
Public Sub RACINE_OUV_DerNum(ByVal nf As Integer, _
                             ByVal a As Double, _
                             ByVal b As Double, _
                             ByVal dx As Double, _
                             ByVal dy As Double, _
                             ByVal xo As Double, _
                             ByVal h As Double, _
                             ByRef X As Double, _
                             ByRef err As Integer, _
                    Optional ByVal M As Integer = 1)
                             
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim d As Double, mh As Double
    
    mh = M * h
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    Y = Fonc(nf, X)
    k = 0: TabX(k) = X
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        d = Fonc(nf, X + h) - Y
        If (d = 0#) Then err = mSINGUL: Exit Sub
        delta = Y * (mh / d)
        X = X - delta: TabX(k) = X
        If prv_Divergence(a, X, b) Then err = mDIVERG: Exit Sub
        Y = Fonc(nf, X)
        TabLamda(k) = Y
        TabDelta(k) = delta
    Loop Until prv_Arret(dx, dy, delta, Y)
    err = k + 1
End Sub
Public Sub RACINE_OUV_Lagrange(ByVal nf As Integer, _
                               ByVal a As Double, _
                               ByVal b As Double, _
                               ByVal dx As Double, _
                               ByVal dy As Double, _
                               ByVal xo As Double, _
                               ByVal h As Double, _
                               ByRef X As Double, _
                               ByRef err As Integer, _
                      Optional ByVal M As Integer = 1)
                            
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim yo As Double
    
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    yo = Fonc(nf, X)
    Y = Fonc(nf, X + h)
    d = yo - Y
    delta = h
    If Abs(Y) > Abs(yo) Then Y = yo Else X = X + h
    k = 0: TabX(k) = X: TabY(k) = Y
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        If (d = 0) Then err = mSINGUL: Exit Sub
        delta = M * Y * (delta / d)
        X = X + delta
        If prv_Divergence(a, X, b) Then err = mDIVERG: Exit Sub
        yo = Y
        Y = Fonc(nf, X)
        d = yo - Y
    Loop Until prv_Arret(dx, dy, delta, Y, k, X)
    err = k + 2
End Sub
Public Sub RACINE_OUV_Muller(ByVal nf As Integer, _
                             ByVal a As Double, _
                             ByVal b As Double, _
                             ByVal dx As Double, _
                             ByVal dy As Double, _
                             ByVal xo As Double, _
                             ByVal h As Double, _
                             ByRef X As Double, _
                             ByRef err As Integer)
                             
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim xn As Double, xp As Double, yo As Double
    Dim yn As Double, yp As Double
    Dim hn As Double, hp As Double
    Dim un As Double, up As Double
    Dim ALPHA As Double, BETA As Double
    Dim GAMMA As Double, w As Double, dum As Double
    
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    xo = X
    yo = Fonc(nf, xo)
    up = h: hp = Fonc(nf, xo + h) - yo
    un = -h: hn = Fonc(nf, xo - h) - yo
    k = 0: TabX(k) = X
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        If un = 0 Or up = 0 Then err = k + 3: Exit Sub
        ALPHA = hp / up - hn / un
        BETA = 0.5 * (up * hn / un - un * hp / up)
        GAMMA = yo * (up - un)
        w = BETA * BETA - ALPHA * GAMMA
        If (w < 0) Then err = mDIVERG: Exit Sub
        w = Sqr(w)
        If BETA < 0 Then w = -w
        w = BETA + w
        delta = -GAMMA / w
        X = xo + delta: TabX(k) = X
        If (b - X) * (X - a) < 0 Then
            delta = -w / ALPHA
            X = xo + delta: TabX(k) = X
            If (b - X) * (X - a) < 0 Then err = mDIVERG: Exit Sub
        End If
        Y = Fonc(nf, X)
        dum = Y - yo
        TabLamda(k) = Y
        TabDelta(k) = delta
        If prv_Arret(dx, dy, delta, Y) Then Exit Do
        If delta > 0 Then
            If delta < up Then
                xo = X
                yo = Y
                un = -delta
                up = up + un
                hn = -dum
                hp = hp + hn
            Else
                xo = xo + up
                yo = yo + hp
                un = -up
                up = delta - up
                hn = -hp
                hp = dum - hn
            End If
        Else
            If delta > un Then
                xo = X
                yo = Y
                up = -delta
                un = un + up
                hp = -dum
                hn = hn + hp
            Else
                xo = xo + un
                yo = yo + hn
                up = -un
                un = delta + up
                hp = -hn
                hn = dum + hp
            End If
        End If
    Loop
    err = k + 3
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
'
'               Autres Méthodes itératives
'
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
Public Sub RACINE_OUV_Multiple(ByVal nf As Integer, _
                               ByVal a As Double, _
                               ByVal b As Double, _
                               ByVal dx As Double, _
                               ByVal dy As Double, _
                               ByVal xo As Double, _
                               ByRef X As Double, _
                               ByRef err As Integer)
                             
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim d As Double, s As Double
    
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    Y = Fonc(nf, X, d, s)
    k = 0: TabX(k) = X
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        delta = d * d - Y * s
        If (delta = 0#) Then err = mSINGUL: Exit Sub
        delta = Y * (d / delta)
        X = X - delta: TabX(k) = X
        If prv_Divergence(a, X, b) Then err = mDIVERG: Exit Sub
        Y = Fonc(nf, X, d, s)
    Loop Until prv_Arret(dx, dy, delta, Y)
    err = k + 1
End Sub
Public Sub RACINE_OUV_Halley(ByVal nf As Integer, _
                             ByVal a As Double, _
                             ByVal b As Double, _
                             ByVal dx As Double, _
                             ByVal dy As Double, _
                             ByVal xo As Double, _
                             ByRef X As Double, _
                             ByRef err As Integer)
                             
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim d As Double, s As Double
    Dim ALPHA As Double
    
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    Y = Fonc(nf, X, d, s)
    k = 0: TabX(k) = X
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        If (d = 0#) Then err = mSINGUL: Exit Sub
        delta = Y / d
        ALPHA = 1 - 0.5 * delta * (s / d)
        If (ALPHA = 0#) Then err = mSINGUL: Exit Sub
        delta = delta / ALPHA
        X = X - delta: TabX(k) = X
        If prv_Divergence(a, X, b) Then err = mDIVERG: Exit Sub
        Y = Fonc(nf, X, d, s)
    Loop Until prv_Arret(dx, dy, delta, Y)
    err = k + 1
End Sub
Public Sub RACINE_OUV_NRmodif(ByVal nf As Integer, _
                              ByVal a As Double, _
                              ByVal b As Double, _
                              ByVal dx As Double, _
                              ByVal dy As Double, _
                              ByVal xo As Double, _
                              ByRef X As Double, _
                              ByRef err As Integer)
                             
    Dim k As Integer
    Dim Y As Double, delta As Double
    Dim d As Double, s As Double
    
    If (a = b) Then err = mNONVAL: Exit Sub
    X = xo
    If prv_Divergence(a, X, b) Then X = 0.5 * (a + b)
    Y = Fonc(nf, X, d, s)
    k = 0: TabX(k) = X
    Do
        k = k + 1
        If (k > mKMAX) Then err = mNONCNV: Exit Sub
        If (d = 0#) Then err = mSINGUL: Exit Sub
        delta = Y / d
        delta = delta * (1 + 0.5 * delta * (s / d))
        X = X - delta: TabX(k) = X
        If prv_Divergence(a, X, b) Then err = mDIVERG: Exit Sub
        Y = Fonc(nf, X, d, s)
    Loop Until prv_Arret(dx, dy, delta, Y)
    err = k + 1
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
'
'          Procédures et fonctions internes
'
' <><><><><><><><><><><><><><><><><><><><><><><><><><>
Private Function prv_Cherche(ByVal nf As Integer, _
                             ByVal a As Double, _
                             ByVal b As Double, _
                             ByRef xn As Double, _
                             ByRef xp As Double, _
                             ByRef yn As Double, _
                             ByRef yp As Double) As Boolean
 
    yp = Fonc(nf, b)
    If (yp > 0#) Then
        xp = b
        xn = a
        yn = Fonc(nf, a)
    Else
        xp = a
        xn = b
        yn = yp
        yp = Fonc(nf, a)
    End If
    prv_Test = (yp * yn >= 0)
 End Function
Private Function prv_Divergence(ByVal a As Double, _
                                ByVal X As Double, _
                                ByVal b As Double) As Boolean

    prv_Divergence = ((X - a) * (b - X) < 0#)
End Function
Private Function prv_Arret(ByVal dx As Double, _
                           ByVal dy As Double, _
                           ByVal delta As Double, _
                           ByVal Y As Double, _
                  Optional ByVal k As Integer = 0, _
                  Optional ByVal X As Double = 0, _
                  Optional ByVal L As Double = 0) As Boolean
 
    prv_Arret = (((Abs(Ecart) < dx) And (Abs(Y) < dy)) Or (Y = 0#))
    TabDelta(k) = delta
    TabY(k) = Y
    TabX(k) = X
    TabLamda(k) = L
End Function
