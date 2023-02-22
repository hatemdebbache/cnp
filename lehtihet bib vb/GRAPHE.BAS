Attribute VB_Name = "GRAPHE"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << GRAPHE >> version 03
'  Requis   :
'  Externes : Fonc(nf,x)
'
'  Descrip. : Graphes 2D de fonctions et données
'
Option Explicit
Sub Init(ByRef frm As Form)
    frm.WindowState = 2
    frm.Show
End Sub
Sub FixeEchelle(ByRef fig As PictureBox, _
                ByVal xmin As Double, _
                ByVal xmax As Double, _
                ByVal ymin As Double, _
                ByVal ymax As Double, _
       Optional ByVal echlog As String)
    
    ' ---- contrôle des données
    Const IDsp = "Echelle"
    Const message = "Vérifier la dimension de la fenêtre"
    If (xmax <= xmin) Or (ymax <= ymin) Then
        Call zprv_Erreur(IDsp, message)
    End If
    ' ---- début du programme
    fig.Cls
    fig.AutoRedraw = True
    fig.ScaleMode = 0
    echlog = UCase(Trim(echlog))
    Select Case echlog
    Case "X", "Y", "XY", "YX"
        If (echlog = "YX") Then echlog = "XY"
        fig.Tag = echlog
        Call zprv_Transf(echlog, xmin, ymin)
        Call zprv_Transf(echlog, xmax, ymax)
    Case Else
        fig.Tag = ""
    End Select
    fig.ScaleLeft = xmin
    fig.ScaleWidth = xmax - xmin
    fig.ScaleTop = ymax
    fig.ScaleHeight = ymin - ymax
End Sub
Sub LireEchelle(ByRef fig As PictureBox, _
                ByRef xmin As Double, _
                ByRef xmax As Double, _
                ByRef ymin As Double, _
                ByRef ymax As Double, _
       Optional ByRef echlog As String)
                     
    xmin = fig.ScaleLeft
    xmax = fig.ScaleWidth + xmin
    ymax = fig.ScaleTop
    ymin = fig.ScaleHeight + ymax
    echlog = fig.Tag
End Sub
Sub TraceAxes(ByRef fig As PictureBox, _
              ByVal xo As Double, _
              ByVal yo As Double, _
              ByVal coulAxes As Integer, _
              ByVal dx As Double, _
              ByVal dy As Double, _
              ByVal coulGrille As Integer, _
     Optional ByVal coulG2 As Integer)
                            
    Dim xmin As Double, xmax As Double, x As Double
    Dim ymin As Double, ymax As Double, y As Double
    Dim coul As Long, coul2 As Long, echlog As String
    Dim k As Integer, z As Double
    Dim Ldix(2 To 10) As Double
    
    ' ---- début du programme
    Call LireEchelle(fig, xmin, xmax, ymin, ymax, echlog)
    If (dx < 0) Then dx = (xmin - xmax) / dx
    If (dy < 0) Then dy = (ymin - ymax) / dy
    If (echlog <> "") Then
        Call zprv_Transf(echlog, xo, yo)
        For k = 2 To 10
            Ldix(k) = zprv_Log(k)
        Next k
    End If
    If (xo < xmin) Then
        xo = xmin
    ElseIf (xo > xmax) Then
        xo = xmax
    Else
    End If
    If (yo < ymin) Then
        yo = ymin
    ElseIf (yo > ymax) Then
        yo = ymax
    Else
    End If
    coul = zprv_coul(coulGrille, fig.ForeColor)
    coul2 = zprv_coul(coulG2, fig.BackColor)
    If (dx > 0) Then
        Select Case echlog
        Case "", "Y"
            For x = xo To xmax Step dx
                fig.Line (x, ymin)-(x, ymax), coul
            Next x
            For x = xo To xmin Step -dx
                fig.Line (x, ymin)-(x, ymax), coul
            Next x
        Case Else
            For x = Fix(xo) To Fix(xmax - 1)
                If (coul2 <> fig.BackColor) Then
                    For k = 2 To 10
                        z = x + Ldix(k)
                        fig.Line (z, ymin)-(z, ymax), coul2
                    Next k
                End If
                fig.Line (x, ymin)-(x, ymax), coul
            Next x
            For x = Fix(xmin) To Fix(xo - 1)
                If (coul2 <> fig.BackColor) Then
                    For k = 2 To 10
                        z = x + Ldix(k)
                        fig.Line (z, ymin)-(z, ymax), coul2
                    Next k
                End If
                fig.Line (x, ymin)-(x, ymax), coul
            Next x
        End Select
    End If
    If (dy > 0) Then
        Select Case echlog
        Case "", "X"
            For y = yo To ymax Step dy
                fig.Line (xmin, y)-(xmax, y), coul
            Next y
            For y = yo To ymin Step -dy
                fig.Line (xmin, y)-(xmax, y), coul
            Next y
        Case Else
            For y = Fix(yo) To Fix(ymax - 1)
                If (coul2 <> fig.BackColor) Then
                    For k = 2 To 10
                        z = y + Ldix(k)
                        fig.Line (xmin, z)-(xmax, z), coul2
                    Next k
                End If
                fig.Line (xmin, y)-(xmax, y), coul
            Next y
            For y = Fix(ymin) To Fix(yo - 1)
                If (coul2 <> fig.BackColor) Then
                    For k = 2 To 10
                        z = y + Ldix(k)
                        fig.Line (xmin, z)-(xmax, z), coul2
                    Next k
                End If
                fig.Line (xmin, y)-(xmax, y), coul
            Next y
        End Select
    End If
    coul = zprv_coul(coulAxes, fig.ForeColor)
    fig.Line (xo, ymin)-(xo, ymax), coul
    fig.Line (xmin, yo)-(xmax, yo), coul
End Sub
Sub TraceFonc(ByRef fig As PictureBox, _
              ByVal nf As Integer, _
              ByVal xmin As Double, _
              ByVal xmax As Double, _
              ByVal Npts As Integer, _
              ByVal relie As Boolean, _
              ByVal coulTrace As Integer, _
     Optional ByVal taille As Integer)
 
    Dim i As Integer
    Dim h As Double, x As Double, y As Double
    Dim coul As Long, siz As Integer
    
    ' ---- contrôle des données
    Const IDsp = "TraceFonc"
    Const message1 = "Npts doit être > 1"
    Const message2 = "xmin doit être différent de xmax"
    If (Npts <= 1) Then Call zprv_Erreur(IDsp, message1)
    If (xmax = xmin) Then Call zprv_Erreur(IDsp, message2)
    ' ---- début du programme
    coul = zprv_coul(coulTrace, fig.ForeColor)
    fig.DrawWidth = zprv_Taille(taille, fig.DrawWidth)
    h = (xmax - xmin) / (Npts - 1)
    x = xmin
    For i = 1 To Npts
        y = Fonc(nf, x)
        Call zprv_Trace(fig, x, y, relie, i - 1, coul)
        x = x + h
    Next i
    fig.DrawWidth = 1
 End Sub
Sub TraceFoncPol(ByRef fig As PictureBox, _
                 ByVal nf As Integer, _
                 ByVal Tmin As Double, _
                 ByVal Tmax As Double, _
                 ByVal Npts As Integer, _
                 ByVal relie As Boolean, _
                 ByVal coulTrace As Integer, _
        Optional ByVal taille As Integer)
 
    Dim i As Integer
    Dim h As Double, T As Double, r As Double
    Dim coul As Long, siz As Integer
    
    ' ---- contrôle des données
    Const IDsp = "TraceFoncPol"
    Const message1 = "Npts doit être > 1"
    Const message2 = "tmin doit être différent de tmax"
    If (Npts <= 1) Then Call zprv_Erreur(IDsp, message1)
    If (Tmax = Tmin) Then Call zprv_Erreur(IDsp, message2)
    ' ---- début du programme
    coul = zprv_coul(coulTrace, fig.ForeColor)
    fig.DrawWidth = zprv_Taille(taille, fig.DrawWidth)
    h = 3.14159265359 * (Tmax - Tmin) / (Npts - 1)
    T = 3.14159265359 * Tmin
    For i = 1 To Npts
        r = Fonc(nf, T)
        Call zprv_Trace(fig, r * Cos(T), r * Sin(T), relie, i - 1, coul)
        T = T + h
    Next i
    fig.DrawWidth = 1
 End Sub
Sub TraceFoncPar(ByRef fig As PictureBox, _
                 ByVal nfx As Integer, _
                 ByVal nfy As Integer, _
                 ByVal Tmin As Double, _
                 ByVal Tmax As Double, _
                 ByVal Npts As Integer, _
                 ByVal relie As Boolean, _
                 ByVal coulTrace As Integer, _
        Optional ByVal taille As Integer)
 
    Dim i As Integer
    Dim h As Double, T As Double, x As Double, y As Double
    Dim coul As Long, siz As Integer
    
    ' ---- contrôle des données
    Const IDsp = "TraceFoncPar"
    Const message1 = "Npts doit être > 1"
    Const message2 = "tmin doit être différent de tmax"
    If (Npts <= 1) Then Call zprv_Erreur(IDsp, message1)
    If (Tmax = Tmin) Then Call zprv_Erreur(IDsp, message2)
    ' ---- début du programme
    coul = zprv_coul(coulTrace, fig.ForeColor)
    fig.DrawWidth = zprv_Taille(taille, fig.DrawWidth)
    h = (Tmax - Tmin) / (Npts - 1)
    T = Tmin
    For i = 1 To Npts
        x = Fonc(nfx, T)
        y = Fonc(nfy, T)
        Call zprv_Trace(fig, x, y, relie, i - 1, coul)
        T = T + h
    Next i
    fig.DrawWidth = 1
 End Sub
Sub TraceVect2(ByRef fig As PictureBox, _
               ByVal Nmin As Integer, _
               ByVal Nmax As Integer, _
               ByRef VecX() As Double, _
               ByRef VecY() As Double, _
               ByVal relie As Boolean, _
               ByVal coulTrace As Integer, _
      Optional ByVal taille As Integer)
                   
    Dim i As Integer
    Dim x As Double, y As Double
    Dim coul As Long
     
    ' ---- début du programme
   coul = zprv_coul(coulTrace, fig.ForeColor)
   fig.DrawWidth = zprv_Taille(taille, fig.DrawWidth)
   For i = Nmin To Nmax
        x = VecX(i)
        y = VecY(i)
        Call zprv_Trace(fig, x, y, relie, i - Nmin, coul)
    Next i
    fig.DrawWidth = 1
End Sub
Sub TraceVect1(ByRef fig As PictureBox, _
               ByVal Nmin As Integer, _
               ByVal Nmax As Integer, _
               ByRef VecY() As Double, _
               ByVal relie As Boolean, _
               ByVal coulTrace As Integer, _
      Optional ByVal taille As Integer)
                   
    Dim i As Integer
    Dim y As Double
    Dim coul As Long
     
    ' ---- début du programme
   coul = zprv_coul(coulTrace, fig.ForeColor)
   fig.DrawWidth = zprv_Taille(taille, fig.DrawWidth)
   For i = Nmin To Nmax
        y = VecY(i)
        Call zprv_Trace(fig, i, y, relie, i - Nmin, coul)
    Next i
    fig.DrawWidth = 1
End Sub
' ##################################
'
'  Procédures et fonctions privées
'
' ##################################
Private Function zprv_coul(ByVal coul As Integer, _
                           ByVal defcoul As Long) As Long

    If (coul >= 0) Then
        zprv_coul = QBColor(coul)
    Else
        zprv_coul = defcoul
    End If
End Function
Private Function zprv_Taille(ByVal taille As Integer, _
                             ByVal defTaille As Integer) As Integer

    If (taille >= 1) Then
        zprv_Taille = taille
    Else
        zprv_Taille = defTaille
    End If
End Function
Private Sub zprv_Trace(ByRef fig As PictureBox, _
                       ByVal x As Double, _
                       ByVal y As Double, _
                       ByVal relie As Boolean, _
                       ByVal k As Integer, _
                       ByVal coul As Long)
                         
    Call zprv_Transf(fig.Tag, x, y)
    If (relie = False) Or (k = 0) Then
        fig.PSet (x, y), coul
    Else
        fig.Line -(x, y), coul
    End If
End Sub
Private Sub zprv_Transf(ByVal echlog As String, _
                        ByRef x As Double, _
                        ByRef y As Double)

    Select Case echlog
    Case "X"
        x = zprv_Log(x)
    Case "Y"
        y = zprv_Log(y)
    Case "XY", "YX"
        x = zprv_Log(x)
        y = zprv_Log(y)
    Case Else
    End Select
End Sub
Private Function zprv_Log(ByVal w As Double) As Double
    
    If (w < 0) Then
        Const message = "L'echelle log ou semi-log exige des données strictement postives"
        Const IDsp = "zprv_Log"
        Call zprv_Erreur(IDsp, message)
    Else
        zprv_Log = 0.4342944819033 * Log(w + 1E-30)
    End If
End Function
Public Sub zprv_Erreur(ByVal IDsp As String, _
                       ByVal message As String)
                            
    Dim p As String
    
    p = "Erreur dans le module GRAPHE " + Chr(13)
    p = p + "Sous-programme : " + IDsp + Chr(13)
    p = p + "Donnée non valide : " + message
    MsgBox p, vbCritical, ""
    Stop
End Sub
                       

